import nanome
from nanome.api.structure import Complex
from nanome.util import async_callback, Logs
from nanome.util.enums import NotificationTypes

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
from dataclasses import dataclass

from .PropertiesHelper import PropertiesHelper


@dataclass
class Entry:
    complex: Complex
    num_frames: int = 0
    ignore_next_update: int = 0

    @property
    def current_frame(self):
        if self.is_conformer:
            return next(self.complex.molecules).current_conformer
        return self.complex.current_frame

    @property
    def current_num_frames(self):
        if self.is_conformer:
            return next(self.complex.molecules).conformer_count
        return len(list(self.complex.molecules))

    @property
    def is_conformer(self):
        return len(list(self.complex.molecules)) == 1


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
        self.server_url = 'data-table-server'
        protocol = 'https' if self.https else 'http'
        try:
            urllib.request.urlopen(f'{protocol}://{self.server_url}')
        except urllib.error.URLError:
            self.server_url = self.url

        self.session = ''.join(random.choices(string.ascii_lowercase, k=4))

        self.selected_entry_index: int = None
        self.selected_indices: list[int] = []
        self.selected_entries: dict[str, Entry] = {}

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
            elif type == 'add-smiles':
                await self.add_smiles(data)
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
        complexes: list[Complex] = await self.request_complex_list()
        items = [{'name': c.full_name, 'index': c.index} for c in complexes]
        await self.ws_send('complexes', items)

        indices = [c.index for c in complexes]
        selected_indices = list(set(self.selected_indices) & set(indices))

        if len(selected_indices) != len(self.selected_indices):
            await self.ws_send('select-complexes', selected_indices)
            await self.select_complexes(selected_indices)

    def on_complex_list_changed(self):
        self.refresh_complexes()

    @async_callback
    async def on_complex_updated(self, complex: Complex):
        Logs.debug('complex updated', complex.index)
        if complex.index not in self.selected_indices:
            return

        entry = self.selected_entries[complex.index]
        entry.complex = complex
        if entry.ignore_next_update:
            entry.ignore_next_update -= 1
            return

        num_frames = entry.current_num_frames
        if num_frames != entry.num_frames:
            entry.num_frames = num_frames
            await self.update_table()
            return

        if complex.index != self.selected_entry_index:
            return

        id = f'{complex.index}-{entry.current_frame}'
        await self.ws_send('select-frame', id)

    def get_entry_and_frame_index(self, id):
        entry_index, frame_index = map(int, id.split('-'))
        entry = self.selected_entries[entry_index]
        return entry, frame_index

    async def select_complexes(self, indices):
        add_indices = list(set(indices) - set(self.selected_indices))
        remove_indices = list(set(self.selected_indices) - set(indices))
        self.selected_indices = indices

        if add_indices:
            complexes: list[Complex] = await self.request_complexes(add_indices)

            for complex in complexes:
                entry = Entry(complex)
                entry.num_frames = entry.current_num_frames
                self.selected_entries[complex.index] = entry
                complex.register_complex_updated_callback(self.on_complex_updated)

        if remove_indices:
            for index in remove_indices:
                del self.selected_entries[index]

            if self.selected_entry_index in remove_indices:
                self.selected_entry_index = None

        await self.update_table()

    async def select_frame(self, id):
        entry, frame = self.get_entry_and_frame_index(id)
        complex = entry.complex
        self.selected_entry_index = complex.index

        if entry.is_conformer:
            next(complex.molecules).set_current_conformer(frame)
        else:
            complex.set_current_frame(frame)

        await self.update_complexes(complex.index)

    async def add_column(self, data, update=True):
        name = data['name']
        values = data['values']
        has_changes = False

        for entry in self.selected_entries.values():
            complex = entry.complex
            if entry.is_conformer:
                molecule = next(complex.molecules)
                for i, associateds in enumerate(molecule.associateds):
                    value = str(values.get(f'{complex.index}-{i}', ''))
                    if value != associateds.get(name):
                        associateds[name] = value
                        has_changes = True
            else:
                for i, molecule in enumerate(complex.molecules):
                    value = str(values.get(f'{complex.index}-{i}', ''))
                    if value != molecule.associated.get(name):
                        molecule.associated[name] = value
                        has_changes = True

        if has_changes and update:
            await self.update_complexes()
            await self.update_table(False)

        return has_changes

    async def add_smiles(self, data):
        smiles = data['smiles']
        hydrogens = data['hydrogens']
        align_to_index = data['align_to']
        align_to = None

        if align_to_index is not None:
            [align_to] = await self.request_complexes([align_to_index])

        existing_complexes = await self.request_complex_list()
        ids = [int(c.name[4:]) for c in existing_complexes if c.name.startswith('VNM-')]
        max_id = max(ids) if ids else 0
        add_complexes = []

        for s in smiles.split('.'):
            complex = self.helper.complex_from_smiles(s, align_to, hydrogens)
            if not complex:
                continue
            complex.name = f'VNM-{max_id + 1:04}'
            add_complexes.append(complex)
            max_id += 1

        if not add_complexes:
            self.send_notification(NotificationTypes.error, "Unable to load invalid structure")
            await self.update_table(False)
            return

        added_complexes = await self.add_to_workspace(add_complexes)
        added_indices = [c.index for c in added_complexes]
        selected_indices = [*self.selected_indices, *added_indices]
        await self.ws_send('select-complexes', selected_indices)
        await self.select_complexes(selected_indices)

    async def delete_frames(self, ids):
        indices = [self.get_entry_and_frame_index(id) for id in ids]
        indices = sorted(indices, key=lambda x: x[1], reverse=True)
        to_update = set()
        to_remove = set()

        for entry, index in indices:
            complex = entry.complex
            entry.num_frames -= 1

            if entry.is_conformer:
                molecule = next(complex.molecules)
                molecule.delete_conformer(index)
            else:
                del complex._molecules[index]

            if entry.num_frames == 0:
                to_remove.add(complex.index)
                to_update.discard(complex.index)
            else:
                to_update.add(complex.index)

        if to_update:
            await self.update_complexes(*to_update)

        if self.selected_entry_index in to_remove:
            self.selected_entry_index = None

        if to_remove:
            remove_complexes = []
            for index in to_remove:
                remove_complexes.append(self.selected_entries[index].complex)
                self.selected_indices.remove(index)
                del self.selected_entries[index]

            await self.remove_from_workspace(remove_complexes)
            await self.refresh_complexes()
            await self.ws_send('select-complexes', self.selected_indices)

        await self.update_table()

    async def split_frames(self, data):
        ids = data['ids']
        single = data['single']
        remove = data['remove']
        name_column = data['name_column']

        indices = [self.get_entry_and_frame_index(id) for id in ids]
        source = None
        source_index = None
        complexes = []

        if single:
            complex = Complex()
            complex.name = indices[0][0].complex.name
            for entry, index in indices:
                source_complex = entry.complex
                if source_complex.index != source_index:
                    source = source_complex.convert_to_frames()
                    source_index = source.index
                complex.add_molecule(source._molecules[index])
            complexes.append(complex)
        else:
            for entry, index in indices:
                source_complex = entry.complex
                if source_complex.index != source_index:
                    source = source_complex.convert_to_frames()
                    source_index = source.index
                complex = Complex()
                molecule = source._molecules[index]
                complex.name = molecule.associated.get(name_column, source_complex.name)
                complex.add_molecule(molecule)
                complexes.append(complex)

        Complex.align_origins(source, *complexes)
        await self.update_structures_deep(complexes)

        if remove:
            await self.delete_frames(ids)

    async def update_frame(self, data):
        entry, index = self.get_entry_and_frame_index(data['id'])
        complex = entry.complex
        del data['id']

        if entry.is_conformer:
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

        complexes: list[Complex] = await self.request_complex_list()
        to_update = []

        for index in indices:
            complex = next((c for c in complexes if c.index == index), None)
            if not complex:
                continue
            entry = self.selected_entries[index]
            entry.complex.position = complex.position
            entry.complex.rotation = complex.rotation
            entry.ignore_next_update += 1
            to_update.append(entry.complex)

        await self.update_structures_deep(to_update)

    async def update_table(self, update_images=True):
        if self.selected_entry_index:
            entry = self.selected_entries[self.selected_entry_index]
            id = f'{entry.complex.index}-{entry.current_frame}'

        complexes = []
        frames = []

        for index, entry in self.selected_entries.items():
            complex = entry.complex.convert_to_frames()
            complex.index = index
            complexes.append(complex)

            for i, mol in enumerate(complex.molecules):
                frames.append({
                    'id': f'{index}-{i}',
                    'frame': f'{complex.full_name} - {i+1}',
                    **mol.associated
                })

        await self.ws_send('frames', frames)

        if self.selected_entry_index:
            await self.ws_send('select-frame', id)

        if update_images:
            for complex in complexes:
                await self.update_images(complex)

    async def calculate_properties(self):
        # property_name -> complex_frame_id -> value
        properties = {}
        has_changes = False
        all_properties = ['SMILES', *self.helper.property_names]

        for entry in self.selected_entries.values():
            complex = entry.complex
            mols = self.helper.complex_to_mols(complex)

            for i, mol in enumerate(mols):
                if mol is None:
                    continue

                id = f'{complex.index}-{i}'
                mol_properties = self.helper.calculate_properties(mol)

                for name in all_properties:
                    if name not in properties:
                        properties[name] = {}
                    properties[name][id] = mol_properties[name]

        for name in all_properties:
            column = {'name': name, 'values': properties.get(name, {})}
            changed = await self.add_column(column, False)
            has_changes = has_changes or changed

        if has_changes:
            await self.update_complexes()
            await self.update_table(False)

    async def update_images(self, complex: Complex):
        try:
            complex.io.to_sdf(self.temp_sdf.name)
            supplier = Chem.SDMolSupplier(self.temp_sdf.name)
        except:
            return

        for i, mol in enumerate(supplier):
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
