import ipywidgets.widgets as widgets
from ipyfilechooser import FileChooser
from ipyfilechooser.utils import normalize_path
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem import SaltRemover
# from rdkit.Chem import Mol
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem
from pathlib import Path
import os
import shutil
from IPython.display import display
import numpy
import re
from contextlib import contextmanager
from functools import partial

import sqlite3

from rdkit import RDLogger



import warnings
warnings.filterwarnings('ignore')
RDLogger.DisableLog('rdApp.info')




## support functions
@contextmanager
def db_connection():
    db_path = Path('//ltvil-freenas5.thermofisher.lt/mzCloudRepo/molfiles.db')
    try:
        con = sqlite3.connect(db_path)
        yield con
    except Exception as e:
        raise e
    finally:
        con.close()

def get_molString(con):
    def inner(compoundId):
        cursor = con.execute("SELECT molstring FROM molstings WHERE compoundId = ?;", (compoundId,))
        data = cursor.fetchone()
        return data
    return inner

    
def sanitize_mol(mol):
    remover = SaltRemover.SaltRemover(defnData=None)
    removed =  remover.StripMol(mol)
    standartizer = rdMolStandardize.Uncharger()
    return standartizer.uncharge(removed)

def string_filter(item, criterion=''):
    if not criterion in item:
        return item

@contextmanager
def navigate(dir):
    current = os.getcwd()
    try:
        os.chdir(dir)
        yield
    except Exception as e:
        print(e)
    finally:
        os.chdir(current)

## business logic

class OwnFileChooser(FileChooser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        target_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'Output'))
        if not os.path.exists(target_path):
            os.mkdir(target_path)
        self._default_path = normalize_path(target_path)
        

class Application:
    def __init__(self) -> None:
        self.data = ''
        self.errors = []
        self.compound_ids = [] #tracking list for compounds in the scope
        self.comp_class_repo = '' # data from repository
        self.raw_files = [] # raw file list
        self.ccomx_files = [] # ccomx file list
        self.paths_to_raws = [] 
        self.comp_class = []
        self.sdf = ''
        self.base_path = ''
        self.comp_list = '' # data frame from repository
        self.options = ['----']
        self.mols = [] #id bound structure data
        self.comp_names = []
        self.mol_path = '//ltvil-freenas5.thermofisher.lt/mzCloudRepo/molfiles.db'

    def get_compound_list(self, data, sheet_name):
        self.comp_list = pd.read_excel(data, sheet_name=sheet_name)

    def get_comp_classes_master(self, data, sheet_name):
        self.comp_class_repo = pd.read_excel(data, sheet_name=sheet_name)

    def get_actual_compound_data_paths(self, plate_name):
        path = Path(self.base_path) / plate_name
        for item in os.listdir(path):
            path_to_item = path / item
            if path_to_item.is_dir():
                cmp_id, *_ = item.split('_')
                self.compound_ids.append(cmp_id)
                self.paths_to_raws.append({
                    'cID': cmp_id,
                    'path': path_to_item,
                })

    def get_raw_files(self):
        with out1:
            print('Collecting raw file info')
        for element in self.paths_to_raws:
            errors = {'raw': False, 'ccomx': False}
            cmp = element.get('cID')
            base_path = element.get('path')
            raws = [item for item in os.listdir(base_path) if item.endswith('.raw')]
            ccomxs = [item for item in os.listdir(base_path) if item.endswith('.ccomx')]
            unmerged_filter = partial(string_filter, criterion='unmerged')
            ccomxs = filter(unmerged_filter, ccomxs)
            if raws:
                for raw in raws:
                    self.raw_files.append({'cID':cmp, 'path': base_path / raw})
            else:
                errors['raw'] = True
                self.errors.append(f'WARNING {cmp}: has no raw files in the folder')
            if ccomxs:
                for ccomx in ccomxs:
                    self.ccomx_files.append({'cID': cmp, 'path': base_path / ccomx})
            else:
                errors['ccomx'] = True
                self.errors.append(f'WARNING {cmp}: has no ccomx files in the folder')
            
            if all(errors.values()):
                self.compound_ids.remove(cmp)
                self.errors.append(f'CRITICAL {cmp}: has no ccomx files in the folder')

                
    def get_compound_names(self):
        with out1:
            print('Collecting compound data')
        for item in self.compound_ids:
            name, *_ = self.comp_list.Name.loc[
                self.comp_list.CompoundID==item
                ].values
            iupac_name, *_ = self.comp_list.IUPACName.loc[
                self.comp_list.CompoundID==item
                ].values
            if isinstance(iupac_name, numpy.float64):
                iupac_name = ''
            if not name:
                app.errors.append(
                    f'CRITICAL {item}: no assigned name, will not be uploaded'
                    )
                app.compound_ids.remove(item)
            else:
                app.comp_names.append({
                'cID': item,
                'name': name,
                'iupac_name': iupac_name, 
            })

    def get_structure_from_string(self, compID):
        mol_string, *_ = self.comp_list.InchiSmiles.loc[self.comp_list.CompoundID==compID].values
        if mol_string:
            if mol_string.startswith('InChI'):
                mol = Chem.MolFromInchi(mol_string)
            else:
                mol = Chem.MolFromSmiles(mol_string)
            return mol
        else:
            app.errors.append(f'INFO {compID}:- No string type molecular data found')
            raise ValueError('No structure data available to compound')

    def get_structure_data(self):
        with out1:
            print('Collecting structural data')
        for item in self.compound_ids:
            mol = ''
            try: 
                with db_connection() as con:
                    mol_grabber = get_molString(con)
                    molstring = mol_grabber(item)
                    mol = Chem.MolFromMolBlock(molstring)
            except OSError as e:
                self.errors.append(f'INFO {item}:- Mol not found, reverting to InChI / smiles')
                try:
                    mol = self.get_structure_from_string(item)
                except ValueError as e:
                    self.compound_ids.remove(item)
                    self.errors.append(f'CRITICAL {item}: is skipped no structure data')
            if mol:
                self.mols.append({
                    'cID': item,
                    'ROMol':sanitize_mol(mol),
                    })
                
    def get_assigned_classes(self):
        with out1:
            print('Collecting class info')
        patt = re.compile(r'\d+')
        for item in self.compound_ids:
            c_cls_obj, *_ = self.comp_list.CompoundClass.loc[
                self.comp_list.CompoundID==item
            ].values
            c_cls = re.findall(patt, str(c_cls_obj))
            if not c_cls:
                c_cls = ['0'] # has to be string
            for cls in c_cls:
                self.comp_class.append(
                    {'cID': item, 'cclass_id': cls}
                )
    def to_dataframes(self):
        obj_list = ['ccomx_files', 'raw_files', 'mols', 'comp_class', 'comp_names']
        for item in obj_list:
            obj = self.__getattribute__(item)
            obj = pd.DataFrame(obj)
            self.__setattr__(item, obj)
    
    def prepare_sdf_compatible_df(self):
        self.sdf = self.comp_names.merge(self.mols, on='cID')

    def prepare_df_file_names_map(self):
        self.raw_name_map = self.comp_names.merge(self.raw_files, on='cID')
        self.raw_name_map['path'] = self.raw_name_map.path.astype(str)
        self.raw_name_map['file'] = self.raw_name_map.path\
            .str.split('\\').str.get(-1)
        self.ccomx_name_map = self.comp_names.merge(self.raw_files, on='cID')
        self.ccomx_name_map['path'] = self.ccomx_files.path.astype(str)
        self.ccomx_name_map['file'] = self.ccomx_name_map.path\
            .str.split('\\').str.get(-1)
    
    def prepare_compound_class_map(self):
        data = self.comp_names.merge(self.comp_class, on='cID')
        data['cclass_id'] = data['cclass_id'].astype(int)
        self.comp_cclass_map = data.merge(self.comp_class_repo, on='cclass_id')
    
    def copy_datafiles(self, filetype):
        current = Path(os.getcwd())
        os.mkdir(filetype)
        destination = current / filetype
        attribute = f'{filetype}_name_map'
        data = self.__getattribute__(attribute)
        for i in range(len(data)):
            filename = data.file.loc[i]
            target = destination / filename
            source = data.path.loc[i]
            print(f'Copying {filename}')
            shutil.copyfile(source, target)

app  = Application()

## ui elements

plate_selector = widgets.Dropdown(options=app.options)
file_chooser = OwnFileChooser('', select_desc='Output File location', icon='fa-floppy-disk')
prepare_data = widgets.Button(description='Prepare data', icon='fa-check')
write_data = widgets.Button(description='Write Files', icon='fa-file')
out1 = widgets.Output()
out2 = widgets.Output()
placeholder_widget = widgets.HTML('<div style="height: 20px;"></div>')
master_list = widgets.FileUpload(description='UploadMasterList')
read_initial_data = widgets.Button(description='ReadData', icon='fa-book')
data_path = widgets.Textarea(placeholder='Base path to data')

# call backs and mechanics

def initial_data_callback(dummy):
    """ registers to read_initial data button """
    if data_path.value:
        app.base_path = data_path.value.strip()
        plate_selector.options = [item for item in os.listdir(data_path.value)]
    else:
        raise ValueError('Base path to data is not specified')
    if master_list.value:
        data = master_list.value
        file, *_ = data.keys()
        data = data[file]['content']
        app.get_compound_list(data, 'masterlist')
        app.get_comp_classes_master(data, 'classes')
    else:
        raise ValueError('No file was uploaded')
    
    with out1:
        print('Initial data read successful')

read_initial_data.on_click(initial_data_callback)



def prepare_data_callback(dummy):
    """ registers with prepare data button """
    with out1:
        print('---Preparing Data---')
    app.get_actual_compound_data_paths(plate_selector.value)
    app.get_raw_files()
    app.get_structure_data()
    app.get_compound_names()
    app.get_assigned_classes()
    app.to_dataframes()
    app.prepare_sdf_compatible_df()
    app.prepare_df_file_names_map()
    app.prepare_compound_class_map()
    with out1:
        print('-Data preparation complete-')
    with out2:
        print('Errors during data preparation:')
        for err in app.errors:
            print(err)


prepare_data.on_click(prepare_data_callback)


def write_data_callback(dummy):
    out_dir = file_chooser.value
    if not Path(out_dir).exists():
        os.mkdir(out_dir)
    with navigate(out_dir):
        with out1:
            print('---Writing SDF---')
        PandasTools.WriteSDF(
            app.sdf, 
            'structures.sdf', 
            idName='cID',
            properties=['cID', 'name'],)
        with out1:    
            print("---Writing CSVs---")
        app.comp_cclass_map[
            ['cID', 'name', 'CompoundClass']
            ].to_csv('compound_classes.csv', sep=';', index=False)
        app.raw_name_map[
            ['cID', 'name', 'file']
            ].to_csv('raw_file_mapping.csv', sep=';', index=False)
        app.ccomx_name_map[
            ['cID', 'name', 'file']
            ].to_csv('ccomx_file_mapping.csv', sep=';', index=False)
        with out1:
            print('---Copying data files---')
            print('--raw files--')
        app.copy_datafiles('raw')
        with out1:
            print('--ccomx files--')
        app.copy_datafiles('ccomx')
        with out1:
            print('<==ALL DONE==>')

write_data.on_click(write_data_callback)

# UI implementation

with out1:
    print('### STATUS MESSAGGES ###')
with out2:
    print('### ERROR MESSAGGES ###')

do_first = widgets.HBox([widgets.VBox([master_list, read_initial_data]), data_path, placeholder_widget])
vbox1 = widgets.VBox([plate_selector, prepare_data, file_chooser, write_data])
vbox2 = widgets.VBox([out1])
vbox3 = widgets.VBox([out2])
hbox = widgets.HBox([vbox1, vbox2, vbox3])
display(do_first, hbox)
