from flask import *
import pandas as pd
import numpy as np
from fileinput import filename
import sqlite3
import os.path

application = Flask(__name__)

application.config['SECRET_KEY'] = 'string'

connection = sqlite3.connect('database.db', check_same_thread=False)
con = connection.cursor()

con.execute('''DROP TABLE IF EXISTS user_inp''')

con.execute('''CREATE TABLE IF NOT EXISTS user_inp(
            natoms INTEGER,
            ncycles INTEGER,
            label1 TEXT,
            label2 TEXT,
            indexstart INTEGER,
            indexend INTEGER,
            filename TEXT,
            natom1 INTEGER
            )''')

def get_db_conn():
    connect = sqlite3.connect('database.db')
    conn = connect.cursor()
    return conn

@application.route("/")
def home():
	return render_template('index.html')

@application.route("/upl", methods = ['GET', 'POST'])	
def uploaded() -> "html":
    if request.method == 'POST':
        f=request.files['file']
        f.save(f.filename)
        n_of_atom_1 = request.form['n_of_atom_1']
        n_of_atom_1 = int(n_of_atom_1)
        get_pre_dists(f.filename, n_of_atom_1)
        with sqlite3.connect('database.db') as con:
            cur = con.cursor()
            load_atoms = cur.execute('''SELECT natoms FROM user_inp''').fetchall()
            n_of_atoms = [p[0] for p in load_atoms][0]
            load_opts = cur.execute('''SELECT ncycles FROM user_inp''').fetchall()
            opt_lines = [x[0] for x in load_opts][0] 
            load_label1 = cur.execute('''SELECT label1 FROM user_inp''').fetchall()
            label1 = [l1[0] for l1 in load_label1]
            load_label2 = cur.execute('''SELECT label2 FROM user_inp''').fetchall()
            label2 = [l2[0] for l2 in load_label2]
            print(type(label1))
        return render_template('upl.html', fname = f.filename, n_of_atom_1=n_of_atom_1,
        n_of_atoms=n_of_atoms, opt_lines=opt_lines, label1=label1, label2=label2)

@application.route("/res", methods = ['GET', 'POST'])
def results() -> "html":
    if request.method == 'POST':
        cycle_n = request.form['cycle_n']
        atom_label_1 = request.form['atom_label_1']
        atom_label_1 = atom_label_1.upper()
        atom_label_2 = request.form['atom_label_2']
        atom_label_2 = atom_label_2.upper()
        res_dist = request.form['res_dist']
        with sqlite3.connect('database.db') as con:
            cur = con.cursor()
            load_filename = cur.execute('''SELECT filename FROM user_inp''').fetchall()
            filename = [a[0] for a in load_filename][0]
            load_natom1 = cur.execute('''SELECT natom1 FROM user_inp''').fetchall()
            n_of_atom_1 = [b[0] for b in load_natom1][0]
        get_dists(filename, n_of_atom_1, cycle_n, atom_label_1, atom_label_2, res_dist)
        return render_template('res.html', dist_list_init=dist_list_init, dist_list=dist_list)

def get_pre_dists(filename, n_of_atom_1):
    with open(filename) as outfile:
        lines = outfile.readlines()
    n_atoms_starter = 'xyz'
    n_atoms_ender = 'end of input'
    for line_index, line in enumerate(lines):
        line = line.lower().strip()
        if n_atoms_starter in line:
            line_index_starter = line_index + 1
        elif n_atoms_ender in line:
            line_index_ender = line_index
    n_of_atoms = line_index_ender - line_index_starter
    opt_lines = 0
    for index, line in enumerate(lines):
	    line = line.lower().strip()
	    if 'geometry optimization cycle' in line:
	        opt_lines += 1
    xyz_init_list_clear = []
    xyz_init_list = lines[line_index_starter:line_index_ender]
    for element in xyz_init_list:
        element = element.replace('|', ' ').replace('>', ' ').replace('\n', ' ')
        xyz_init_list_clear.append(element.strip().split())
    xyz_init_list_1 = xyz_init_list_clear[:n_of_atom_1]
    xyz_init_list_2 = xyz_init_list_clear[n_of_atom_1:]
    xyz_init_df_1 = pd.DataFrame(xyz_init_list_1, columns=['atom index', 'atom', 'x','y','z'])
    xyz_init_df_2 = pd.DataFrame(xyz_init_list_2, columns=['atom index', 'atom', 'x','y','z'])
    xyz_init_df_1.loc[:, 'atom'] = xyz_init_df_1.loc[:, 'atom'].astype(str)
    xyz_init_df_2.loc[:, 'atom'] = xyz_init_df_2.loc[:, 'atom'].astype(str)
    atom_label_1_pos = xyz_init_df_1.atom.unique()
    atom_label_2_pos = xyz_init_df_2.atom.unique()  
    label1str = f"{','.join(atom_label_1_pos)}"
    label2str = f"{','.join(atom_label_2_pos)}"
    with sqlite3.connect('database.db') as con:
        cur = con.cursor()
        con.execute('''INSERT INTO user_inp VALUES(?,?,?,?,?,?,?,?)''', 
        (n_of_atoms, opt_lines, label1str, label2str, line_index_starter, line_index_ender, filename, n_of_atom_1))
    
def get_dists(filename, n_of_atom_1, cycle_n, atom_label_1, atom_label_2, res_dist):
    load_atoms = con.execute('''SELECT natoms FROM user_inp''').fetchall()
    n_of_atoms = [p[0] for p in load_atoms][0]
    load_opts = con.execute('''SELECT ncycles FROM user_inp''').fetchall()
    opts = [x[0] for x in load_opts][0] 
    load_label1 = con.execute('''SELECT label1 FROM user_inp''').fetchall()
    label1 = [l1[0] for l1 in load_label1]
    load_label2 = con.execute('''SELECT label2 FROM user_inp''').fetchall()
    label2 = [l2[0] for l2 in load_label2]
    load_start_index = con.execute('''SELECT indexstart FROM user_inp''').fetchall()
    load_end_index = con.execute('''SELECT indexend FROM user_inp''').fetchall()
    start_index = int([y[0] for y in load_start_index][0])
    end_index = int([z[0] for z in load_end_index][0])
    with open(filename) as outfile:
        lines = outfile.readlines()  
    global dist_list_init
    global dist_list
    dist_list_init, dist_list = [],[]
    cycle_n = int(cycle_n)
    res_dist = float(res_dist)
    if cycle_n > opts:
        print("Entered number is out of range")
    elif cycle_n > 99:
        cycle_string = f"cycle {cycle_n}"
    elif cycle_n > 9:
        cycle_string = f"cycle  {cycle_n}"
    else:
        cycle_string = f"cycle   {cycle_n}"
    xyz_init_list_clear = []
    xyz_init_list = lines[start_index:end_index]
    for element in xyz_init_list:
        element = element.replace('|', ' ').replace('>', ' ').replace('\n', ' ')
        xyz_init_list_clear.append(element.strip().split())
    xyz_init_list_1 = xyz_init_list_clear[:n_of_atom_1]
    xyz_init_list_2 = xyz_init_list_clear[n_of_atom_1:]
    xyz_init_df_1 = pd.DataFrame(xyz_init_list_1, columns=['atom index', 'atom', 'x','y','z'])
    xyz_init_df_2 = pd.DataFrame(xyz_init_list_2, columns=['atom index', 'atom', 'x','y','z'])
    xyz_init_df_1.loc[:, ['x','y','z']] = xyz_init_df_1.loc[:, ['x','y','z']].astype(float)
    xyz_init_df_1.loc[:, 'atom'] = xyz_init_df_1.loc[:, 'atom'].astype(str)
    xyz_init_df_2.loc[:, ['x','y','z']] = xyz_init_df_2.loc[:, ['x','y','z']].astype(float)
    xyz_init_df_2.loc[:, 'atom'] = xyz_init_df_2.loc[:, 'atom'].astype(str)
    init_chs_1 = xyz_init_df_1['atom'].isin([atom_label_1])
    init_chs_2 = xyz_init_df_2['atom'].isin([atom_label_2])
    xyz_init_df_1_chs = xyz_init_df_1[init_chs_1]
    xyz_init_df_2_chs = xyz_init_df_2[init_chs_2]
    for i in range(len(xyz_init_df_1_chs.index)):
	    c = np.array([xyz_init_df_1_chs.iloc[i,2],xyz_init_df_1_chs.iloc[i,3],xyz_init_df_1_chs.iloc[i,4]])
	    for k in range(len(xyz_init_df_2_chs)):
	        d = np.array([xyz_init_df_2_chs.iloc[k,2],xyz_init_df_2_chs.iloc[k,3],xyz_init_df_2_chs.iloc[k,4]])
	        dist_init = np.linalg.norm(c-d)
	        if dist_init < res_dist:
                 dist_cycle_init = f"Initial Geometry: Distance between {atom_label_1} {xyz_init_df_1_chs.iloc[i,0]} and {atom_label_2} {xyz_init_df_2_chs.iloc[k,0]} atom: {dist_init} nm"
                 dist_list_init.append(dist_cycle_init)
    for index, line in enumerate(lines):
	    line = line.lower().strip()
	    if cycle_string in line:
	        s_index = index +5
	        e_index = s_index + n_of_atoms
	        xyz_list_cycle= lines[s_index:e_index]
    xyz_list_cycle_clear = []
    for elem in xyz_list_cycle:		
        elem = elem.replace('|', ' ').replace('>', ' ').replace('\n', ' ')
        xyz_list_cycle_clear.append(elem.strip().split())
    xyz_list_cycle1 = xyz_list_cycle_clear[:n_of_atom_1]
    xyz_list_cycle2 = xyz_list_cycle_clear[n_of_atom_1:]
    xyz_cycle_df_1 = pd.DataFrame(xyz_list_cycle1, columns=['atom', 'x','y','z'])
    xyz_cycle_df_2 = pd.DataFrame(xyz_list_cycle2, columns=['atom', 'x','y','z'])
    xyz_cycle_df_1['atom index'] = xyz_init_df_1['atom index']
    xyz_cycle_df_2['atom index'] = xyz_init_df_2['atom index']
    xyz_cycle_df_1 = xyz_cycle_df_1[['atom index', 'atom', 'x', 'y', 'z']]
    xyz_cycle_df_2 = xyz_cycle_df_2[['atom index', 'atom', 'x', 'y', 'z']]
    xyz_cycle_df_1.loc[:, ['x','y','z']] = xyz_cycle_df_1.loc[:, ['x','y','z']].astype(float)
    xyz_cycle_df_2.loc[:, ['x','y','z']] = xyz_cycle_df_2.loc[:, ['x','y','z']].astype(float)
    xyz_cycle_df_1.loc[:, 'atom'] = xyz_cycle_df_1.loc[:, 'atom'].astype(str)
    xyz_cycle_df_2.loc[:, 'atom'] = xyz_cycle_df_2.loc[:, 'atom'].astype(str)
    cycle_chs_1 = xyz_cycle_df_1['atom'].isin([atom_label_1])
    cycle_chs_2 = xyz_cycle_df_2['atom'].isin([atom_label_2])
    xyz_cycle_df_1_chs = xyz_cycle_df_1[cycle_chs_1]
    xyz_cycle_df_2_chs = xyz_cycle_df_2[cycle_chs_2]
    for i in range(len(xyz_cycle_df_1_chs.index)):
	    a = np.array([xyz_cycle_df_1_chs.iloc[i,2],xyz_cycle_df_1_chs.iloc[i,3],xyz_cycle_df_1_chs.iloc[i,4]])
	    for k in range(len(xyz_cycle_df_2_chs)):
	        b = np.array([xyz_cycle_df_2_chs.iloc[k,2],xyz_cycle_df_2_chs.iloc[k,3],xyz_cycle_df_2_chs.iloc[k,4]])
	        dist = np.linalg.norm(a-b)
	        if dist < res_dist:
                 dist_cycle = f"Geometry cycle {cycle_n}: Distance between {atom_label_1} {xyz_cycle_df_1_chs.iloc[i,0]} and {atom_label_2} {xyz_cycle_df_2_chs.iloc[k,0]} atom: {dist} nm"
                 dist_list.append(dist_cycle)    
    con.close()
    if os.path.isfile('results.txt'):
        os.remove('results.txt')
    make_txt(cycle_n)
    return dist_list_init, dist_list

def make_txt(cycle_n): 
    with open('results.txt', 'a') as results:
        for elem in dist_list_init:
            results.write(f'\n {elem}')
        for elem in dist_list:
            results.write(f'\n {elem}')
        results.close()

if __name__ == "__main__":
	application.run(debug = True)