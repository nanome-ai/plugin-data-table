import nanome
from nanome.util import async_callback, Logs

from cairosvg import svg2png
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

import argparse
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

class WebTable(nanome.AsyncPluginInstance):
    @async_callback
    async def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)

        self.server_url = self.custom_data[0]
        self.session = '-'.join(''.join(random.choices(string.ascii_lowercase, k=3)) for _ in range(3))
        self.ws = await websockets.connect(f'ws://web-table-server/ws')
        self.selected_complex = None
        self.selected_frame = None
        self.ignore_next_update = False

        await self.ws_send('host', self.session)
        self.on_run()

        while True:
            m = await self.ws.recv()
            msg = json.loads(m)
            data = msg.get('data')

            if msg['type'] == 'join':
                self.update_complexes()
            elif msg['type'] == 'select-complex':
                await self.select_complex(data)
            elif msg['type'] == 'select-frame':
                await self.select_frame(data)

    def on_run(self):
        self.open_url(f'http://{self.server_url}/{self.session}')

    @async_callback
    async def on_stop(self):
        await self.ws.close()

    async def ws_send(self, type, data):
        msg = json.dumps({'type': type, 'data': data})
        await self.ws.send(msg)

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
        if complex.index != self.selected_complex.index:
            return
        if self.ignore_next_update:
            self.ignore_next_update = False
            return

        frame = self.get_selected_frame(complex)
        await self.ws_send('select-frame', frame)

    async def select_complex(self, index):
        complexes = await self.request_complexes([index])
        if len(complexes) == 0:
            return

        complex = complexes[0]
        self.selected_is_conformer = len(list(complex.molecules)) == 1

        self.selected_complex = complex
        complex.register_complex_updated_callback(self.on_complex_updated)
        frame = self.get_selected_frame(complex)

        data = []
        complex = complex.convert_to_frames()
        complex.index = complexes[0].index

        for i, mol in enumerate(complex.molecules):
            data.append({'index': i, **mol.associated})

        await self.ws_send('data', data)
        await self.ws_send('select-frame', frame)
        await self.generate_images(complex)

    async def select_frame(self, index):
        if self.selected_complex is None:
            return

        self.ignore_next_update = True
        if self.selected_is_conformer:
            next(self.selected_complex.molecules).set_current_conformer(index)
        else:
            self.selected_complex.set_current_frame(index)
        await self.update_structures_deep([self.selected_complex])

    def get_selected_frame(self, complex):
        if self.selected_is_conformer:
            return next(complex.molecules).current_conformer
        return complex.current_frame

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
