from nanome.api.structure import Complex
from nanome.util import Logs

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
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

API_CACHE_TIME = 30
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
            Property('AR', '%d', mDesc.CalcNumAromaticRings),
            Property('InChI', '%s', Chem.MolToInchi),
            Property('InChIKey', '%s', Chem.MolToInchiKey),
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
    def property_names(self):
        return [p.name for p in self.properties]

    def fetch_property(self, endpoint, prop, mol):
        name = endpoint['name']
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        cache_id = f'{name}:{smiles}'
        cache = self.api_cache.get(cache_id)

        if cache is None or datetime.now() > cache['expiration']:
            url = endpoint['url']
            method = endpoint['method']
            data = endpoint['data']
            headers = endpoint.get('headers', {})

            try:
                if data == 'smiles' and method == 'GET':
                    url = url.replace(':smiles', quote(smiles))
                    result = requests.get(url, headers=headers).json()

                elif data == 'smiles' and method == 'POST':
                    payload = endpoint['payload']
                    if type(payload) != str:
                        payload = json.dumps(payload)
                    payload = payload.replace(':smiles', smiles)
                    if headers.get('Content-Type') != 'application/json':
                        payload = json.loads(payload)
                    result = requests.post(url, headers=headers, data=payload).json()

                elif data == 'sdf' and method == 'POST':
                    Chem.SDWriter(self.temp_sdf.name).write(mol)
                    files = {'file': open(self.temp_sdf.name, 'rb')}
                    result = requests.post(url, headers=headers, files=files).json()

                else:
                    Logs.error(f'Unsupported request type: {method} {data} on {name}')
                    return None

            except:
                Logs.error(f'Failed to fetch {name}')
                return None

            data = {}
            for item, info in endpoint['properties'].items():
                value = result
                data[item] = None

                try:
                    for path in info['path'].replace('][', '].[').split('.'):
                        match = re.search(r'([^\[]+)?(?:\[(?:(\d+)|(.+?))\])?', path)
                        (key, index, query) = match.groups()

                        if key is not None:
                            value = value.get(key)
                        if index is not None:
                            value = value[int(index)]
                        if query is not None:
                            q_key, q_value = query.split('=', 1)
                            for v in value:
                                if v.get(q_key) == q_value:
                                    value = v
                                    break

                    if type(value) not in [str, int, float, bool]:
                        raise TypeError

                    data[item] = value

                except:
                    Logs.error(f'Invalid path for {item} on {name}')

            cache_time = timedelta(seconds=endpoint.get('cache_time', API_CACHE_TIME))
            cache = { 'expiration': datetime.now() + cache_time, 'data': data }
            self.api_cache[cache_id] = cache

        return cache['data'][prop]

    def calculate_properties(self, mol):
        Chem.AssignStereochemistryFrom3D(mol)
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        cache = self.smiles_to_property_cache.get(smiles)

        if cache:
            return cache

        properties = { 'SMILES': smiles }
        for name, fmt, fn in self.properties:
            value = fn(mol)
            value = 'ERR' if  value is None else fmt % value
            properties[name] = value

        self.smiles_to_property_cache[smiles] = properties
        return properties

    def complex_from_smiles(self, smiles, align_to_complex=None, hydrogens=True):
        mol = Chem.MolFromSmiles(smiles)
        if not smiles or mol is None:
            return None

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)

        if not hydrogens:
            mol = Chem.RemoveHs(mol)

        if align_to_complex is not None:
            align_to_mol = next(self.complex_to_mols(align_to_complex))
            align_to_mol.UpdatePropertyCache()

            # remove hydrogens for alignment
            mol = Chem.RemoveHs(mol)
            align_to_mol = Chem.RemoveHs(align_to_mol)

            atomMap = None
            try:
                _, _, atomMap = AllChem.GetBestAlignmentTransform(mol, align_to_mol)
            except:
                try:
                    # sometimes changing order of molecules helps
                    _, _, atomMap = AllChem.GetBestAlignmentTransform(align_to_mol, mol)
                    atomMap = [(b, a) for a, b in atomMap]
                except:
                    Logs.error(f'Failed to align {Chem.MolToSmiles(mol)} to {Chem.MolToSmiles(align_to_mol)}')

            if atomMap:
                AllChem.AlignMol(mol, align_to_mol, atomMap=atomMap)

            if hydrogens:
                mol = Chem.AddHs(mol)

        with Chem.SDWriter(self.temp_sdf.name) as w:
            w.SetForceV3000(True)
            w.write(mol)

        complex = Complex.io.from_sdf(path=self.temp_sdf.name)
        if align_to_complex is not None:
            complex.position = align_to_complex.position
            complex.rotation = align_to_complex.rotation

        return complex

    def complex_to_mols(self, complex: Complex):
        complex.io.to_sdf(self.temp_sdf.name)
        supplier = Chem.SDMolSupplier(self.temp_sdf.name)
        return supplier

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
