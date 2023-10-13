import ipywidgets.widgets as widgets
from ipyfilechooser import FileChooser
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem import SaltRemover
from rdkit.Chem import Mol
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem
from pathlib import Path
import os
import shutil
from IPython.display import display
from typing import List, Tuple
#import requests
#from requests.auth import HTTPBasicAuth
import re


def debug(func):
    def wrapper(*args, **kwargs):
        rv = func(*args, **kwargs)
        print(f'{func.__name__} ran')
        print('---------')
        print(rv)
        print('---------')
        return rv
    return wrapper


## constansts

# def download_file():
#     local_filename = 'unified_data.xlsx'
#     # NOTE the stream=True parameter below
#     with requests.get(
#         'https://thermofisher.sharepoint.com/:x:/s/PTVG/EUqeF7ld3bpBpQ4mCopMP9ABUyyxSjobz2kqrI_4jK3SFg?e=YIJTOU / user',
#         auth = HTTPBasicAuth('ltvil.ptvg', 'tmo#VATlab@00001'),stream=True) as r:
#         r.raise_for_status()
#         with open(local_filename, 'wb') as f:
#             for chunk in r.iter_content(chunk_size=8192): 
#                 # If you have chunk encoded response uncomment if
#                 # and set chunk_size parameter to None.
#                 #if chunk: 
#                 f.write(chunk)

# download_file()
base_path = '//ltvil-freenas5.thermofisher.lt/PTVG_Data/DATA/mzCloud/Data/Production'

compound_list = 'unified_data.xlsx'
main = {'data': '',
        'errors': [],
        'cclass_list':'',
        'raw_files':'',
        'sdf': ''}

## UI matters

def run_drop_down(base_path=base_path):
    options = [item for item in os.listdir(base_path)]
    options.insert(0, '---')
    return widgets.Dropdown(options=options)

plate_selector = run_drop_down()
file_chooser = FileChooser('', select_desc='Output File location', icon='fa-floppy-disk')
submit = widgets.Button(description='Prepare data', icon='fa-check')
write_data = widgets.Button(description='Write Files', icon='fa-file')
out1 = widgets.Output()
out2 = widgets.Output()

## business logic
# @debug
def get_compound_table(filters: List[str]=None) ->pd.DataFrame:
    compound_data = pd.read_excel(compound_list)
    if filters:
        compound_data = compound_data.loc[compound_data.ID.isin(filters)].copy()
    return compound_data.reset_index()

# @debug
def get_actual_compound_data_paths(path: Path) ->Tuple[List[str], List[str]]:
    paths_to_raws = []
    compound_ids = []
    for item in os.listdir(path):
        path_to_item = path / item
        if path_to_item.is_dir():
            paths_to_raws.append(path_to_item)
            cmp_id, *_ = item.split('_')
            # print(cmp_id)
            compound_ids.append(cmp_id)
    # i = 0
    # while True:
    #     try:
    #         print(compound_ids[i], '\t', paths_to_raws[i])
    #         i += 1
    #     except IndexError as e:
    #         print(e)
    #         break
    # print(len(set(paths_to_raws)), len(set(compound_ids)))
    return paths_to_raws, compound_ids

# @debug
# def get_actual_compound_ids(paths: List[str]) -> List[str]:
#     to_return = []
#     for path in paths:
#         cmp_id, *_ = path.split('_')
#         to_return.append(cmp_id)
#     return to_return

#@debug
def get_mols() ->List[Mol]:
    global main
    mols=[]
    for item in main['data'].InChiSmiles.values:
        if 'InChI' in item: 
            mol = Chem.MolFromInchi(item)
        else:
            mol = Chem.MolFromSmiles(item)
        if isinstance(mol, Mol):
            mols.append(mol)
        else:
            print("shit happened", item)
    return mols

#@debug
def remove_salt(mol: Mol) -> str:
    remover = SaltRemover.SaltRemover(defnData=None)
    removed =  remover.StripMol(mol)
    standartizer = rdMolStandardize.Uncharger()
    return standartizer.uncharge(removed)

#@debug
def mol_to_inchi(mol: Mol) -> str:
    return Chem.MolToInchi(mol)

#@debug
def mol_to_smiles(mol: Mol) -> str:
    return Chem.MolToSmiles(mol)

# @debug
def make_data_frame(target_path: Path) -> pd.DataFrame:
    act_cmp_paths, act_cmp_ids = get_actual_compound_data_paths(target_path)
    # act_cmp_ids = get_actual_compound_ids(act_cmp_paths)

    data_frame = get_compound_table(filters=act_cmp_ids)
    data_frame['to_drop'] = "False"
    data_frame['base_path'] = act_cmp_paths
    return data_frame[['ID', 'Name', 'CompClass', 'InChiSmiles', 'base_path', 'to_drop']]


def add_meta_data() -> pd.DataFrame:
    global main
    data = []
    search_pat = re.compile(r'\d+')
    main['data']['CompClass'] = main['data'].CompClass.astype(str)
    for i in range(len(main['data'])):
        cmp_id = main['data'].ID.loc[i]
        name = main['data'].Name.loc[i]
        inchi = main['data'].InChi.loc[i]
        classes = main['data'].CompClass.loc[i]
        classes = re.findall(search_pat, classes)
        for cls in classes:
            data.append({
                'ID': cmp_id,
                'Name':name,
                'InChi':inchi,
                'class_id': int(cls)
            })
    data_frame = pd.DataFrame(data)
    db = pd.read_excel('classes.xlsx', sheet_name='Sheet1')
    return data_frame.merge(right=db, how='inner', on='class_id')

def collect_raw_files():
    global main
    errors = []
    raw_file_list = []
    for i in range(len(main['data'])):
        id = main['data'].ID.loc[i]
        inchi = main['data'].InChi.loc[i]
        path = main['data'].base_path.loc[i]
        name = main['data'].Name.loc[i]
        raws = [file for file in os.listdir(Path(path)) if file.endswith('.ccomx')]
        if len(raws) == 0:
            errors.append({"ID": id, 'Error': "No Raw file"})
            main['data'].to_drop.at[i] = True
        else:
            for raw in raws:
                raw_file_list.append(
                    {
                        'Name': name,
                        'ID': id,
                        'InChI': inchi,
                        'file': raw,
                        'raw_path': path / raw,
                    }
                )
    return pd.DataFrame(raw_file_list), pd.DataFrame(errors)

# def sanity_check():
#     global main
#     raw_files = s']
#     errors = main['errors']
#     for cmp in errors:
#         if cmp in raw_files['Compound_ID']:
#             raw_files.drop(raw_files['Compound_ID']==cmp, axis=0, inplace=True)
#     raw_files.reset_index(inplace=True)
#     main['raw_files'] = raw_files 

def copy_ccomx():
    global main
    os.mkdir('ccomx')
    current = os.getcwd()
    base = Path(current) / 'ccomx'
    raw_files = main['raw_files']
    for i in range(len(raw_files)):
        file = raw_files.file[i]
        with out1:
            print(f'{file}')
        
        source = Path(raw_files.raw_path.loc[i])
        dest = base / file
        shutil.copyfile(source, dest)

def make_sdf_compatible():
    df = main['data'][['ID', 'Name', 'Smiles', 'InChi']].copy()
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol='Smiles', molCol='ROMol')
    return df

def write_sdf(sdf_data, out_file):
    PandasTools.WriteSDF(sdf_data, out_file, idName='ID', properties=['ID', 'Name', 'InChi', 'Smiles'])
    print("file written")

### mechanics:

def write_files(unused):
    global main
    current_dir = os.getcwd()
    if not file_chooser.value:
        with out2:
            print('Choose the directory')
    else:
        os.mkdir(file_chooser.value)
        os.chdir(file_chooser.value)
        main['cclass_list'].to_csv('class_list.csv', sep=',', index=False)
        main['raw_files'].to_csv('raw_file_list.csv', sep=',', index=False)
        with out1:
            print ('--- Writing SDF file---')
        write_sdf(main['data'], out_file='compound_list.sdf')
        with out1:
            print('--- Copying raw files---')
        copy_ccomx()
    with out1:
        print('ALL DONE')

 
def prepare_data(unused):
    global main
    path = Path(base_path) / plate_selector.value
    with out1:
        print('---Running data collection---')
    main['data'] = make_data_frame(path)
    mols = get_mols()
    mols = list(map(remove_salt, mols))

    updated_inchi = list(map(mol_to_inchi, mols))
    smiless = list(map(mol_to_smiles, mols))
    main['data']['InChi'] = ''
    main['data']['Smiles'] = ''
    for i, items in enumerate(zip(updated_inchi, smiless)):
        inchi, smiles = items
        main['data'].InChi.at[i] = inchi
        main['data'].Smiles.at[i] = smiles
    main['data'].drop('InChiSmiles', axis=1, inplace=True)
    PandasTools.AddMoleculeColumnToFrame(main['data'], smilesCol='Smiles')
    with out1:
        print('---Data collection finnished---')
    with out1:
        print('---Collecting Raw files---')
    main['raw_files'], errors = collect_raw_files()
    main['errors'].append(errors)
    with out2:
        print(main['errors'][0])
    with out1:
        print('---Raw file collection finished---')
    with out1:
        print('---Getting metadata---')
    main['cclass_list'] = add_meta_data()
    with out1:
        print('---Finished Getting metadata---')
    sdf = make_sdf_compatible()
    # with out1:
    #     print('---Performing clean-up---')
    # sanity_check()
    # with out1:
    #     print('---Data ready to export---')
    

submit.on_click(prepare_data)
write_data.on_click(write_files)




### display

with out1:
    print('######STATUS MESSAGGES######')
with out2:
    print('######ERROR MESSAGGES######')

vbox1 = widgets.VBox([plate_selector, submit, file_chooser, write_data])
vbox2 = widgets.VBox([out1])
vbox3 = widgets.VBox([out2])
hbox = widgets.HBox([vbox1, vbox2, vbox3])
display(hbox)
