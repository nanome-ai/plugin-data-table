import nanome
from nanome.api.structure import Complex
from nanome.util import async_callback, Logs

from cairosvg import svg2png
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

import argparse
import asyncio
import base64
import json
import os
import random
import string
import tempfile
import websockets

# mol 2d image drawing options
Draw.DrawingOptions.atomLabelFontSize = 40
Draw.DrawingOptions.dotsPerAngstrom = 100
Draw.DrawingOptions.bondLineWidth = 8

IS_DOCKER = os.path.exists('/.dockerenv')

class WebTable(nanome.AsyncPluginInstance):
    @async_callback
    async def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)

        self.url = self.custom_data[0]
        self.server_url = 'web-table-server' if IS_DOCKER else self.url
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
        ws_url = f'ws://{self.server_url}/ws'
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
            elif type == 'delete-frames':
                await self.delete_frames(data)
            elif type == 'reorder-frames':
                await self.reorder_frames(data)
            elif type == 'select-complex':
                await self.select_complex(data)
            elif type == 'select-frame':
                await self.select_frame(data)
            elif type == 'split-frames':
                await self.split_frames(data)

            Logs.debug('done', type)

    def on_run(self):
        self.open_url(f'{self.url}/{self.session}')

    @async_callback
    async def on_stop(self):
        await self.ws.close()

    @async_callback
    async def update_complexes(self):
        complexes = await self.request_complex_list()
        items = [{'name': c.full_name, 'index': c.index} for c in complexes]
        await self.ws_send('complexes', items)

    def on_complex_added(self):
        self.update_complexes()

    def on_complex_removed(self):
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

    async def reorder_frames(self, indices):
        if self.selected_is_conformer:
            molecule = next(self.selected_complex.molecules)
            num_frames = molecule.conformer_count
            for i, index in enumerate(indices):
                molecule.copy_conformer(index, num_frames + i)
            for _ in range(num_frames):
                molecule.delete_conformer(0)
        else:
            molecules = list(self.selected_complex.molecules)
            self.selected_complex._molecules = [molecules[i] for i in indices]
            for molecule in self.selected_complex.molecules:
                molecule.index = -1

        await self.update_complex()
        await self.update_table()

    async def split_frames(self, data):
        indices = data['indices']
        single = data['single']
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
        await self.delete_frames(indices)

    async def update_complex(self):
        self.ignore_next_update += 1
        await self.update_structures_deep([self.selected_complex])

    async def update_table(self):
        frame = self.get_selected_frame(self.selected_complex)

        complex = self.selected_complex.convert_to_frames()
        complex.index = self.selected_complex.index
        self.selected_num_frames = len(list(complex.molecules))
        data = []

        for i, mol in enumerate(complex.molecules):
            data.append({'index': i, **mol.associated})

        await self.ws_send('data', data)
        await self.ws_send('select-frame', frame)
        await self.generate_images(complex)

    async def generate_images(self, complex, frame=None):
        try:
            complex.io.to_sdf(self.temp_sdf.name)
            supplier = Chem.SDMolSupplier(self.temp_sdf.name)
        except:
            return

        if len(supplier) == 0:
            return

        frames = [frame] if frame is not None else range(len(supplier))

        for i in frames:
            mol = supplier[i]
            if mol is None:
                continue

            Chem.AssignStereochemistryFrom3D(mol)
            AllChem.Compute2DCoords(mol)
            mol = Draw.rdMolDraw2D.PrepareMolForDrawing(mol)

            width, height = 256, 192
            drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
            drawer.drawOptions().additionalAtomLabelPadding = 0.3
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            svg = svg.replace('stroke-linecap:butt', 'stroke-linecap:round')

            id = f'{complex.index}-{i}'
            path = os.path.join(self.temp_dir.name, f'{id}.png')
            svg2png(bytestring=svg, write_to=path, output_width=width, output_height=height)

            await self.send_image(id)

    async def send_image(self, id):
        path = os.path.join(self.temp_dir.name, f'{id}.png')
        with open(path, 'rb') as f:
            data = base64.b64encode(f.read()).decode('utf-8')
            await self.ws_send('image', {'id': id, 'data': data})

def main():
    parser = argparse.ArgumentParser(description='Parse arguments for Web Table plugin')
    parser.add_argument('-u', '--url', dest='url', type=str, help='URL of the web server', required=True)
    args, _ = parser.parse_known_args()

    plugin = nanome.Plugin('Web Table', 'A Nanome plugin to demonstrate the ability to integrate with in-VR web browser', 'other', False)
    plugin.set_plugin_class(WebTable)
    plugin.set_custom_data(args.url)
    plugin.run()


if __name__ == '__main__':
    main()
