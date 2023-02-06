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

        self.selected_complex = None
        self.selected_frame = None
        self.selected_is_conformer = False
        self.selected_num_frames = 0
        self.ignore_next_update = 0

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
                self.update_complexes()
            elif type == 'add-column':
                await self.add_column(data)
            elif type == 'calculate-properties':
                await self.calculate_properties()
            elif type == 'delete-frames':
                await self.delete_frames(data)
            elif type == 'select-complex':
                await self.select_complex(data)
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
    async def update_complexes(self):
        complexes = await self.request_complex_list()
        items = [{'name': c.full_name, 'index': c.index} for c in complexes]
        await self.ws_send('complexes', items)

    def on_complex_list_changed(self):
        self.update_complexes()

    @async_callback
    async def on_complex_updated(self, complex):
        Logs.debug('complex updated', complex.index)
        if complex.index != self.selected_complex.index:
            return

        self.selected_complex = complex
        if self.ignore_next_update:
            self.ignore_next_update -= 1
            return

        if self.selected_is_conformer:
            num_frames = next(complex.molecules).conformer_count
        else:
            num_frames = len(list(complex.molecules))

        if num_frames != self.selected_num_frames:
            await self.select_complex(complex.index)
        else:
            frame = self.get_selected_frame(complex)
            await self.ws_send('select-frame', frame)

    def get_selected_frame(self, complex):
        if self.selected_is_conformer:
            return next(complex.molecules).current_conformer
        return complex.current_frame

    async def select_complex(self, index):
        complexes = await self.request_complexes([index])
        if not complexes:
            return

        complex = complexes[0]
        self.selected_complex = complex
        self.selected_is_conformer = len(list(complex.molecules)) == 1
        complex.register_complex_updated_callback(self.on_complex_updated)

        await self.update_table()

    async def select_frame(self, index):
        if self.selected_is_conformer:
            next(self.selected_complex.molecules).set_current_conformer(index)
        else:
            self.selected_complex.set_current_frame(index)

        await self.update_complex()

    async def add_column(self, data, update=True):
        name = data['name']
        values = data['values']

        if self.selected_is_conformer:
            molecule = next(self.selected_complex.molecules)
            for i, associateds in enumerate(molecule.associateds):
                associateds[name] = str(values[i])
        else:
            for i, molecule in enumerate(self.selected_complex.molecules):
                molecule.associated[name] = str(values[i])

        if update:
            await self.update_complex()
            await self.update_table(False)

    async def delete_frames(self, indices):
        indices = sorted(indices, reverse=True)
        if self.selected_is_conformer:
            molecule = next(self.selected_complex.molecules)
            for index in indices:
                molecule.delete_conformer(index)
        else:
            for index in indices:
                del self.selected_complex._molecules[index]

        await self.update_complex()
        await self.update_table()

    async def split_frames(self, data):
        indices = data['indices']
        single = data['single']
        remove = data['remove']
        name_column = data['name_column']
        source = self.selected_complex.convert_to_frames()
        complexes = []

        if single:
            complex = Complex()
            complex.name = source.name
            for index in indices:
                complex.add_molecule(source._molecules[index])
            complexes.append(complex)
        else:
            for index in indices:
                complex = Complex()
                molecule = source._molecules[index]
                complex.name = molecule.associated[name_column]
                complex.add_molecule(molecule)
                complexes.append(complex)

        Complex.align_origins(source, *complexes)
        await self.update_structures_deep(complexes)

        if remove:
            await self.delete_frames(indices)

    async def update_frame(self, data):
        index = data['index']
        del data['index']

        if self.selected_is_conformer:
            molecule = next(self.selected_complex.molecules)
            for name, value in data.items():
                molecule.associateds[index][name] = value
        else:
            molecule = self.selected_complex._molecules[index]
            molecule.associated.update(data)

        await self.update_complex()

    async def update_complex(self):
        complexes = await self.request_complex_list()
        complex = next(c for c in complexes if c.index == self.selected_complex.index)
        self.selected_complex.position = complex.position
        self.selected_complex.rotation = complex.rotation
        self.ignore_next_update += 1
        await self.update_structures_deep([self.selected_complex])

    async def update_table(self, update_images=True):
        frame = self.get_selected_frame(self.selected_complex)

        complex = self.selected_complex.convert_to_frames()
        complex.index = self.selected_complex.index
        self.selected_num_frames = len(list(complex.molecules))
        frames = []

        for i, mol in enumerate(complex.molecules):
            frames.append({'index': i, **mol.associated})

        await self.ws_send('frames', frames)
        await self.ws_send('select-frame', frame)

        if update_images:
            await self.update_images(complex)

    async def calculate_properties(self):
        try:
            self.selected_complex.io.to_sdf(self.temp_sdf.name)
            supplier = Chem.SDMolSupplier(self.temp_sdf.name)
        except:
            return

        properties = {}
        for name in self.helper.property_names:
            properties[name] = [''] * len(supplier)

        for i, mol in enumerate(supplier):
            if mol is None:
                continue

            mol_properties = self.helper.calculate_properties(mol)
            for name in self.helper.property_names:
                properties[name][i] = mol_properties[name]

        for name in self.helper.property_names:
            await self.add_column({'name': name, 'values': properties[name]}, False)

        await self.update_complex()
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
