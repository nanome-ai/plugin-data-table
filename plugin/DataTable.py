import nanome
from nanome.api.structure import Complex
from nanome.util import async_callback, Logs

from rdkit import Chem

import argparse
import asyncio
import json
import os
import random
import string
import tempfile
import urllib.request
import websockets

from .PropertiesHelper import PropertiesHelper


class DataTable(nanome.AsyncPluginInstance):

    @async_callback
    async def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)

        self.helper = PropertiesHelper(self.temp_dir)

        self.url, self.https = self.custom_data
        self.ws = None

        # If 'data-table-server' not available on local docker network,
        # use external url as server_url
        priority_server_host = 'data-table-server'
        protocol = 'https' if self.https else 'http'
        try:
            urllib.request.urlopen(f'{protocol}://{priority_server_host}')
        except urllib.error.URLError:
            self.server_url = self.url
        else:
            self.server_url = priority_server_host

        self.session = ''.join(random.choices(string.ascii_lowercase, k=4))

        self.selected_complex_index = None
        self.selected_indices = []
        self.selected_complexes = {}
        self.complex_is_conformer = {}
        self.complex_num_frames = {}
        self.ignore_next_update = {}

        await self.ws_connect()
        self.on_run()
        self.ws_loop()

    async def ws_connect(self):
        protocol = 'wss' if self.https else 'ws'
        ws_url = f'{protocol}://{self.server_url}/ws'
        Logs.debug(f'connecting to {ws_url}')

        while True:
            try:
                self.ws = await websockets.connect(ws_url)
                break
            except:
                await asyncio.sleep(1)

        Logs.debug(f'connected to {ws_url}')
        await self.ws_send('host', self.session)

    async def ws_send(self, type, data):
        msg = json.dumps({'type': type, 'data': data})
        await self.ws.send(msg)

    @async_callback
    async def ws_loop(self):
        while True:
            try:
                m = await asyncio.wait_for(self.ws.recv(), timeout=0.1)
            except asyncio.TimeoutError:
                continue
            except websockets.exceptions.ConnectionClosedError:
                await self.ws_connect()
                continue

            msg = json.loads(m)
            type = msg.get('type')
            data = msg.get('data')

            Logs.debug('recv', type, data)

            if type == 'join':
                self.refresh_complexes()
            elif type == 'add-column':
                await self.add_column(data)
            elif type == 'calculate-properties':
                await self.calculate_properties()
            elif type == 'delete-frames':
                await self.delete_frames(data)
            elif type == 'select-complexes':
                await self.select_complexes(data)
            elif type == 'select-frame':
                await self.select_frame(data)
            elif type == 'split-frames':
                await self.split_frames(data)
            elif type == 'update-frame':
                await self.update_frame(data)

            Logs.debug('done', type)

    def on_run(self):
        protocol = 'https' if self.https else 'http'
        self.open_url(f'{protocol}://{self.url}/{self.session}')

    @async_callback
    async def on_stop(self):
        del self.helper
        if self.ws:
            await self.ws.close()

    @async_callback
    async def refresh_complexes(self):
        complexes = await self.request_complex_list()
        items = [{'name': c.full_name, 'index': c.index} for c in complexes]
        await self.ws_send('complexes', items)

    def on_complex_list_changed(self):
        self.refresh_complexes()

    @async_callback
    async def on_complex_updated(self, complex):
        Logs.debug('complex updated', complex.index)
        if complex.index not in self.selected_indices:
            return

        self.selected_complexes[complex.index] = complex
        if self.ignore_next_update[complex.index]:
            self.ignore_next_update[complex.index] -= 1
            return

        num_frames = self.get_num_frames(complex)
        if num_frames != self.complex_num_frames[complex.index]:
            self.complex_num_frames[complex.index] = num_frames
            await self.update_table()
            return

        if complex.index != self.selected_complex_index:
            return

        frame = self.get_selected_frame(complex)
        id = f'{complex.index}-{frame}'
        await self.ws_send('select-frame', id)

    def get_complex_frame(self, id):
        complex_index, frame_index = map(int, id.split('-'))
        return self.selected_complexes[complex_index], frame_index

    def get_num_frames(self, complex):
        if self.complex_is_conformer[complex.index]:
            return next(complex.molecules).conformer_count
        return len(list(complex.molecules))

    def get_selected_frame(self, complex):
        if self.complex_is_conformer[complex.index]:
            return next(complex.molecules).current_conformer
        return complex.current_frame

    async def select_complexes(self, indices):
        add_indices = list(set(indices) - set(self.selected_indices))
        remove_indices = list(set(self.selected_indices) - set(indices))
        self.selected_indices = indices

        if add_indices:
            complexes = await self.request_complexes(add_indices)

            for complex in complexes:
                self.selected_complexes[complex.index] = complex
                self.complex_is_conformer[complex.index] = len(list(complex.molecules)) == 1
                self.complex_num_frames[complex.index] = self.get_num_frames(complex)
                self.ignore_next_update[complex.index] = 0
                complex.register_complex_updated_callback(self.on_complex_updated)

        if remove_indices:
            for index in remove_indices:
                del self.selected_complexes[index]
                del self.complex_is_conformer[index]
                del self.complex_num_frames[index]
                del self.ignore_next_update[index]

            if self.selected_complex_index in remove_indices:
                self.selected_complex_index = None

        await self.update_table()

    async def select_frame(self, id):
        complex, frame = self.get_complex_frame(id)
        self.selected_complex_index = complex.index

        if self.complex_is_conformer[complex.index]:
            next(complex.molecules).set_current_conformer(frame)
        else:
            complex.set_current_frame(frame)

        await self.update_complexes(complex.index)

    async def add_column(self, data, update=True):
        name = data['name']
        values = data['values']

        for complex in self.selected_complexes.values():
            if self.complex_is_conformer[complex.index]:
                molecule = next(complex.molecules)
                for i, associateds in enumerate(molecule.associateds):
                    id = f'{complex.index}-{i}'
                    associateds[name] = str(values[id])
            else:
                for i, molecule in enumerate(complex.molecules):
                    id = f'{complex.index}-{i}'
                    molecule.associated[name] = str(values[id])

        if update:
            await self.update_complexes()
            await self.update_table(False)

    async def delete_frames(self, ids):
        indices = [self.get_complex_frame(id) for id in ids]
        indices = sorted(indices, key=lambda x: x[1], reverse=True)
        to_update = set()
        to_remove = set()

        for complex, index in indices:
            self.complex_num_frames[complex.index] -= 1

            if self.complex_is_conformer[complex.index]:
                molecule = next(complex.molecules)
                molecule.delete_conformer(index)
            else:
                del complex._molecules[index]

            if self.complex_num_frames[complex.index] == 0:
                to_remove.add(complex.index)
                to_update.discard(complex.index)
            else:
                to_update.add(complex.index)

        if to_update:
            await self.update_complexes(*to_update)

        if to_remove:
            remove_complexes = []
            for index in to_remove:
                remove_complexes.append(self.selected_complexes[index])
                self.selected_indices.remove(index)
                del self.selected_complexes[index]
                del self.complex_is_conformer[index]
                del self.complex_num_frames[index]
                del self.ignore_next_update[index]

            await self.remove_from_workspace(remove_complexes)
            await self.refresh_complexes()
            await self.ws_send('select-complexes', self.selected_indices)

        await self.update_table()

    async def split_frames(self, data):
        ids = data['ids']
        single = data['single']
        remove = data['remove']
        name_column = data['name_column']

        indices = [self.get_complex_frame(id) for id in ids]
        source = None
        source_index = None
        complexes = []

        if single:
            complex = Complex()
            complex.name = source.name
            for source_complex, index in indices:
                if source_complex.index != source_index:
                    source = source_complex.convert_to_frames()
                    source_index = source.index
                complex.add_molecule(source._molecules[index])
            complexes.append(complex)
        else:
            for source_complex, index in indices:
                if source_complex.index != source_index:
                    source = source_complex.convert_to_frames()
                    source_index = source.index
                complex = Complex()
                molecule = source._molecules[index]
                complex.name = molecule.associated[name_column]
                complex.add_molecule(molecule)
                complexes.append(complex)

        Complex.align_origins(source, *complexes)
        await self.update_structures_deep(complexes)

        if remove:
            await self.delete_frames(ids)

    async def update_frame(self, data):
        complex, index = self.get_complex_frame(data['id'])
        del data['id']

        if self.complex_is_conformer[complex.index]:
            molecule = next(complex.molecules)
            for name, value in data.items():
                molecule.associateds[index][name] = value
        else:
            molecule = complex._molecules[index]
            molecule.associated.update(data)

        await self.update_complexes(complex.index)

    async def update_complexes(self, *indices):
        if not indices:
            indices = self.selected_indices

        complexes = await self.request_complex_list()
        to_update = []

        for index in indices:
            complex = next(c for c in complexes if c.index == index)
            self.selected_complexes[index].position = complex.position
            self.selected_complexes[index].rotation = complex.rotation
            self.ignore_next_update[index] += 1
            to_update.append(self.selected_complexes[index])

        await self.update_structures_deep(to_update)

    async def update_table(self, update_images=True):
        if self.selected_complex_index:
            complex = self.selected_complexes[self.selected_complex_index]
            frame = self.get_selected_frame(complex)
            id = f'{complex.index}-{frame}'

        complexes = []
        frames = []

        for index, complex in self.selected_complexes.items():
            complex = complex.convert_to_frames()
            complex.index = index
            complexes.append(complex)

            for i, mol in enumerate(complex.molecules):
                frames.append({
                    'id': f'{index}-{i}',
                    'frame': f'{complex.full_name} - {i+1}',
                    **mol.associated
                })

        await self.ws_send('frames', frames)

        if self.selected_complex_index:
            await self.ws_send('select-frame', id)

        if update_images:
            for complex in complexes:
                await self.update_images(complex)

    async def calculate_properties(self):
        # property_name -> complex_frame_id -> value
        properties = {}

        for complex in self.selected_complexes.values():
            try:
                complex.io.to_sdf(self.temp_sdf.name)
                supplier = Chem.SDMolSupplier(self.temp_sdf.name)
            except:
                return

            for i, mol in enumerate(supplier):
                if mol is None:
                    continue

                id = f'{complex.index}-{i}'
                mol_properties = self.helper.calculate_properties(mol)

                for name in self.helper.property_names:
                    if name not in properties:
                        properties[name] = {}
                    properties[name][id] = mol_properties[name]

        for name in self.helper.property_names:
            await self.add_column({'name': name, 'values': properties[name]}, False)

        await self.update_complexes()
        await self.update_table(False)

    async def update_images(self, complex):
        try:
            complex.io.to_sdf(self.temp_sdf.name)
            supplier = Chem.SDMolSupplier(self.temp_sdf.name)
        except:
            return

        if len(supplier) == 0:
            return

        for i in range(len(supplier)):
            mol = supplier[i]
            if mol is None:
                continue

            id = f'{complex.index}-{i}'
            data = self.helper.render_image(mol)
            await self.ws_send('image', {'id': id, 'data': data})


def main():
    parser = argparse.ArgumentParser(description='Parse arguments for Data Table plugin')
    parser.add_argument('--https', dest='https', action='store_true', help='Enable HTTPS on the Data Table Web UI')
    parser.add_argument('-u', '--url', dest='url', type=str, help='URL of the web server')
    parser.add_argument('-w', '--web-port', dest='web_port', type=int, help='Custom port for connecting to Data Table Web UI.')
    args, _ = parser.parse_known_args()

    https = args.https or os.environ.get('SERVER_HTTPS', False)
    port = args.web_port or os.environ.get('SERVER_PORT', None)
    url = args.url or os.environ.get('SERVER_URL', '')

    if not url:
        raise Exception('--url flag or SERVER_URL environment variable must be set')

    if port:
        url = f'{url}:{port}'

    plugin = nanome.Plugin('Data Table', 'A Nanome plugin to view multi-frame structure metadata in a table', 'Analysis', False)
    plugin.set_plugin_class(DataTable)
    plugin.set_custom_data(url, https)
    plugin.run()


if __name__ == '__main__':
    main()
