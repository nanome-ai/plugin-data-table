from nanome.util import Logs

from rdkit import Chem
from rdkit.Chem import Draw
import rdkit.Chem.Descriptors as Desc
import rdkit.Chem.rdMolDescriptors as mDesc

import base64
import json
import os
import re
import requests
import tempfile
from cairosvg import svg2png
from collections import namedtuple
from datetime import datetime, timedelta
from functools import partial
from urllib.parse import quote

# mol 2d image drawing options
Draw.DrawingOptions.atomLabelFontSize = 40
Draw.DrawingOptions.dotsPerAngstrom = 100
Draw.DrawingOptions.bondLineWidth = 8

API_CACHE_TIME = timedelta(seconds=1)
API_SETTINGS = os.path.join(os.path.dirname(__file__), '..', 'config.json')

Property = namedtuple('Property', ['name', 'format', 'fn'])

class PropertiesHelper:
    def __init__(self, temp_dir):
        self.temp_dir = temp_dir
        self.temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)

        self.api_cache = {}
        self.smiles_to_property_cache = {}
        self.smiles_to_image_cache = {}

        self.properties = [
            Property('MW', '%.3f', Desc.MolWt),
            Property('logP', '%.3f', lambda mol: mDesc.CalcCrippenDescriptors(mol)[0]),
            Property('TPSA', '%.3f', mDesc.CalcTPSA),
            Property('HBA', '%d', mDesc.CalcNumHBA),
            Property('HBD', '%d', mDesc.CalcNumHBD),
            Property('RB', '%d', mDesc.CalcNumRotatableBonds),
            Property('AR', '%d', mDesc.CalcNumAromaticRings)
        ]

        if not os.path.exists(API_SETTINGS):
            return

        with open(API_SETTINGS, 'r') as f:
            self.api = json.load(f)

        if self.api.get('overwrite'):
            self.properties = []

        # validate config
        try:
            required_endpoint_keys = ['url', 'method', 'data']
            required_property_keys = ['format', 'path']

            for endpoint in self.api.get('endpoints'):
                if not endpoint.get('name'):
                    raise Exception('Invalid config: missing endpoint name')

                test = [endpoint[k] for k in required_endpoint_keys]
                if endpoint['method'] == 'POST' and endpoint['data'] == 'smiles':
                    test = endpoint['payload']
                for prop, info in endpoint.get('properties').items():
                    test = [info[k] for k in required_property_keys]

        except KeyError as key:
            raise Exception(f'Invalid config: missing {key} on {endpoint["name"]}')

        except TypeError:
            raise Exception(f'Invalid config: array where object should be')

        # register properties
        for endpoint in self.api.get('endpoints'):
            for prop, info in endpoint['properties'].items():
                fn = partial(self.fetch_property, endpoint, prop)
                p = Property(prop, info['format'], fn)
                self.properties.append(p)

    @property
    def num_props(self):
        return len(self.properties)

    @property
    def property_names(self):
        return [p.name for p in self.properties]

    def fetch_property(self, endpoint, prop, mol):
        name = endpoint['name']
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        cache_id = f'{name}:{smiles}'
        cache = self.api_cache.get(cache_id)

        if cache is None or datetime.now() - cache['time'] > API_CACHE_TIME:
            url = endpoint['url']
            method = endpoint['method']
            data = endpoint['data']

            try:
                if data == 'smiles' and method == 'GET':
                    url = url.replace(':smiles', quote(smiles))
                    json = requests.get(url).json()

                elif data == 'smiles' and method == 'POST':
                    payload = endpoint['payload'].replace(':smiles', smiles)
                    headers = {'Content-Type': 'application/json'}
                    json = requests.post(url, headers=headers, data=payload).json()

                elif data == 'sdf' and method == 'POST':
                    Chem.SDWriter(self.temp_sdf.name).write(mol)
                    files = {'file': open(self.temp_sdf.name, 'rb')}
                    json = requests.post(url, files=files).json()

                else:
                    Logs.error(f'Unsupported request type: {method} {data} on {name}')
                    return None

            except:
                Logs.error(f'Failed to fetch {name}')
                return None

            data = {}
            for item, info in endpoint['properties'].items():
                value = json
                data[item] = None

                try:
                    for path in info['path'].replace('][', '].[').split('.'):
                        match = re.search(r'([^\[]+)?(?:\[(\d+)\])?', path)
                        (key, index) = match.groups()

                        if key is not None:
                            value = value.get(key)
                        if index is not None:
                            value = value[int(index)]

                    if type(value) not in [str, int, float, bool]:
                        raise TypeError

                    data[item] = value

                except:
                    Logs.error(f'Invalid path for {item} on {name}')

            cache = { 'time': datetime.now(), 'data': data }
            self.api_cache[cache_id] = cache

        return cache['data'][prop]

    def calculate_properties(self, mol):
        Chem.AssignStereochemistryFrom3D(mol)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        cache = self.smiles_to_property_cache.get(smiles)

        if cache:
            return cache

        properties = {}
        for name, fmt, fn in self.properties:
            value = fn(mol)
            value = 'ERR' if  value is None else fmt % value
            properties[name] = value

        self.smiles_to_property_cache[smiles] = properties
        return properties

    def render_image(self, mol):
        Chem.AssignStereochemistryFrom3D(mol)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        cache = self.smiles_to_image_cache.get(smiles)

        if cache:
            return cache

        Chem.rdCoordGen.AddCoords(mol)
        mol = Draw.rdMolDraw2D.PrepareMolForDrawing(mol)

        width, height = 256, 192
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
        options = drawer.drawOptions()
        Draw.rdMolDraw2D.SetDarkMode(options)
        options.additionalAtomLabelPadding = 0.2
        options.clearBackground = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svg = svg.replace('stroke-linecap:butt', 'stroke-linecap:round')

        png = tempfile.NamedTemporaryFile(delete=False, suffix='.png', dir=self.temp_dir.name)
        svg2png(bytestring=svg, write_to=png.name, output_width=width*2, output_height=height*2)

        with open(png.name, 'rb') as f:
            data = base64.b64encode(f.read()).decode('utf-8')
            self.smiles_to_image_cache[smiles] = data

        return data
