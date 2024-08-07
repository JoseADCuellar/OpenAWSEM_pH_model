#####Definitivo####
import numpy as np
import pandas as pd
import random
import math

def convert_sequence_to_three_letter(seq):
    # Diccionario de mapeo de una letra a tres letras
    one_to_three = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
        'C': 'CYS', 'E': 'GLU', 'Q': 'GLN', 'G': 'GLY',
        'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
        'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
        'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
    }
    # Convertir la secuencia de una letra a tres letras
    seq_list = [one_to_three[aa] for aa in seq]
    return seq_list

def elegir_residuo_mc(lista_residuos):
    residuo_mc, charge = random.choice(lista_residuos)
    return residuo_mc, charge

def aob_residuo_mc(residuo_mc, type_dict, seq_list):
    nombre_residuo_mc = seq_list[residuo_mc]     
    tipo_residuo_mc = 1 if type_dict[nombre_residuo_mc] == 'B' else -1
    return tipo_residuo_mc, nombre_residuo_mc

def calcular_distancia(punto1, punto2):
    return np.sqrt((punto1[0] - punto2[0])**2 + (punto1[1] - punto2[1])**2 + (punto1[2] - punto2[2])**2)

def calcular_distancias_residuo(residuo_mc, cb_positions_df):
    # Obtener la posición del residuo de interés
    posicion_interes = cb_positions_df[cb_positions_df['Residue_Number'] == residuo_mc][['X', 'Y', 'Z']].values[0]
    
    # Calcular las distancias a todos los demás residuos
    vecinos = []
    distancias = []
    for index, row in cb_positions_df.iterrows():
        if row['Residue_Number'] != residuo_mc:
            posicion_actual = [row['X'], row['Y'], row['Z']]
            distancia = calcular_distancia(posicion_interes, posicion_actual)
            vecinos.append(row['Residue_Number'])  # Cambio de nombre a 'vecinos'
            distancias.append(distancia)
    
    return vecinos, distancias

def obtener_nombres_residuos(vecinos, seq_list):
    vecinos = [int(vecino) for vecino in vecinos]
    # Obtener los nombres de los residuos correspondientes a los índices en 'vecinos'
    residue_names = [seq_list[residuo] for residuo in vecinos]
    return residue_names

def generate_Vcharge(residue_sequence_numbers, lista_residuos):
    Vcharge = []
    for residuo in residue_sequence_numbers:
        charge = next((carga for r, carga in lista_residuos if r == residuo), 0)
        Vcharge.append(charge)
    return Vcharge

def clasificar_vecinos_aob(residue_names, type_dict):
    Acidos = [-1 if type_dict[res] == 'A' else 0 for res in residue_names]
    Basicos = [1 if type_dict[res] == 'B' else 0 for res in residue_names]
    Glicinas = [1 if res == 'GLY' else 0 for res in residue_names]
    return Acidos, Basicos, Glicinas

def conteo_polares_no_polares(residue_names, distancias, type_dict):
    # Parametros Polares y No_polares, radio de la esfera para medir
    R_max = 5/10  # Umbral para distancia de polares
    R_max_no_polares = 7/10  # Umbral para distancia de no polares
    tau = 0.1*100  # Valor para tau en la fórmula
    
    Polares = []
    No_polares = []
    for res, distancia in zip(residue_names, distancias):
        if type_dict[res] in ['P', 'B', 'A']:
            valor_polares = 1 if distancia <= R_max else math.exp(-tau * (distancia - R_max)**2)
            Polares.append(valor_polares)
        else:
            Polares.append(0)
    
        if type_dict[res] == 'N':
            valor_no_polares = 1 if distancia <= R_max_no_polares else math.exp(-tau * (distancia - R_max_no_polares)**2)
            No_polares.append(valor_no_polares)
        else:
            No_polares.append(0)

    return Polares, No_polares

def calcular_self(Resname, aob, Num_polares, Num_no_polares):
    if Resname == 'ASP':
        Bp = 0.72426
        Bnp = 5.57075
    elif Resname == 'GLU':
        Bp = 0.380763
        Bnp = 9.44846
    elif Resname == 'LYS' or Resname == 'ARG':
        Bp = 0.0110311
        Bnp = 5.65882
    elif Resname == 'HIS':
        Bp = 0.602896
        Bnp = 9.36507
    elif Resname == 'CYS':
        Bp = 0.037 * 0.001987 * 300
        Bnp = 85.91 * 0.001987 * 300
    elif Resname == 'TYR':
        Bp = 0.0
        Bnp = 0.0

    Npmax = 3.9
    Nnpmax = 17.78
    alfaP = 0.416683
    alfanp = 0.049
    
    # Calculate Up and Unp
    Up = math.exp(-alfaP * (Num_polares - Npmax)**2) if Num_polares <= Npmax else 1
    Unp = math.exp(-alfanp * (Num_no_polares - Nnpmax)**2) if Num_no_polares <= Nnpmax else 1
    
    # Calculate the term
    term_self = aob  * (Bnp * Unp - Bp * Up)
    
    return term_self


def calcular_elec(data, charge):
    K_elec = 2.43232/10
    L = 10.0/10.0
    E_elect = []
    
    for index, row in data.iterrows():
        vecinos_carga = row['Vcharge']
        distancia = row['Distancia']

        if vecinos_carga == -1:  # Vecino ácido
            Calcular = charge * ((-1 / distancia) * math.exp(-distancia / L))
            E_elect.append(Calcular)
        elif vecinos_carga == 1:  # Vecino básico
            Calcular = charge * ((1 / distancia) * math.exp(-distancia / L))
            E_elect.append(Calcular)
        else:  # Vecino sin carga
            E_elect.append(0)

    data['E_elect'] = E_elect
    # Sumar la energia Electrostatica
    term_elec = data['E_elect'].sum() * K_elec
    return term_elec

def suma_polares_no_polares(data):
    polares = data['Polar'].sum()  # Sumar vecinos polares
    no_polares = data['No_Polar'].sum()  # Sumar vecinos no polares
    return no_polares, polares

def charge_flip(aob, old_charge):
    new_charge = 0.0
    if aob == -1.0:  # Ácido
        if old_charge == -1.0:
            new_charge = 0.0
        elif old_charge == 0.0:
            new_charge = -1.0
    elif aob == 1:  # Base
        if old_charge == 1.0:
            new_charge = 0.0
        elif old_charge == 0.0:
            new_charge = 1.0
    return new_charge
def estatico(matriz_posiciones, charge_residues, seq):
    type_dict = {
        'ALA': 'N', 'ARG': 'B', 'ASN': 'P', 'ASP': 'A', 'CYS': 'A', 'GLU': 'A', 'GLN': 'P', 'GLY': 'G', 'HIS': 'B',
        'ILE': 'N', 'LEU': 'N', 'LYS': 'B', 'MET': 'N', 'PHE': 'N', 'PRO': 'N', 'SER': 'P', 'THR': 'P', 'TRP': 'N',
        'TYR': 'A', 'VAL': 'N'
    } 

    aminoacidos_cargados = charge_residues
    posiciones = matriz_posiciones
    seq_list = convert_sequence_to_three_letter(seq)  # Convierte la secuencia a lista de códigos de tres letras

    # Elegir residuo al azar
    residuo_mc, charge = elegir_residuo_mc(aminoacidos_cargados)
    # Ver aob de residuo
    tipo_residuo_mc, nombre_residuo_mc = aob_residuo_mc(residuo_mc, type_dict, seq_list)
    # Ver los vecinos y sus distancias
    vecinos, distancias = calcular_distancias_residuo(residuo_mc, posiciones)
    # Obtener nombres de los vecinos
    residue_names = obtener_nombres_residuos(vecinos, seq_list)
    # Generar Vcharge
    Vcharge = generate_Vcharge(vecinos, charge_residues)
    # Clasificar vecinos
    Acidos, Basicos, Glicinas = clasificar_vecinos_aob(residue_names, type_dict)
    # Contar polares y no polares
    Polares, No_polares = conteo_polares_no_polares(residue_names, distancias, type_dict)

    # Armar un DataFrame con la información
    data = pd.DataFrame({
        'Residuo_mc': residuo_mc,
        'nombre_residuo_mc': nombre_residuo_mc,
        'Carga': charge,
        'Vecinos': vecinos,
        'Vecinos_name': residue_names,
        'Vcharge': Vcharge,
        'Acido': Acidos,
        'Base': Basicos,
        'Polar': Polares,
        'No_Polar': No_polares,
        'Glicina': Glicinas,
        'Distancia': distancias
    })

    # Calcular energías
    no_polares, polares = suma_polares_no_polares(data)
    term_self = calcular_self(nombre_residuo_mc, tipo_residuo_mc, polares, no_polares)
    term_elec = calcular_elec(data, charge)

    # Resultados
    Resultado = {
        'Residuo_mc': residuo_mc,
        'nombre_residuo_mc': nombre_residuo_mc,
        'carga': charge,
        'aob': tipo_residuo_mc,
        'Num_polares': polares,
        'Num_no_polares': no_polares,
        'Term_self': term_self,
        'Term_elec': term_elec
    }
    return pd.DataFrame(Resultado, index=[0]), data, charge_residues
    
def Montecarlo(ph, lista_residuos, cb_positions_df, seq, enlaces_matrix):
    # Parámetros
    kb = 0.001987
    T = 300
    
    lista_cambios = []
#    with open('charge.txt', 'r') as file:
#        lista_residuos = []
#        for line in file.readlines():
#            residuo, carga = line.split()
#            carga = float(carga)
#            if carga != 0.0:
#                lista_residuos.append((int(residuo), carga)) 
    
    old, data, _ = estatico(cb_positions_df, lista_residuos, seq)
    
    elec = old['Term_elec'].iloc[0]
    self = old['Term_self'].iloc[0]
    Num_polares = old['Num_polares'].iloc[0]
    Num_no_polares = old['Num_no_polares'].iloc[0]
    old_carga = old['carga'].iloc[0]
    Resname = old['nombre_residuo_mc'].iloc[0]
    residuo_mc = old['Residuo_mc'].iloc[0]
    aob = old['aob'].iloc[0]

    new_carga = charge_flip(aob, old_carga)
    delta_carga = new_carga - old_carga
        
    if Resname == 'ASP':
        Pka_ref = 4.0
    elif Resname == 'GLU':
        Pka_ref = 4.5
    elif Resname == 'LYS':
        Pka_ref = 10.6
    elif Resname == 'ARG':
        Pka_ref = 12.0
    elif Resname == 'HIS':
        Pka_ref = 6.4
    elif Resname == 'CYS':
    	Pka_ref = 8.3
    elif Resname == 'TYR':
    	Pka_ref = 11.0

    # Calcular nuevos términos de energía
    new_elec = calcular_elec(data, new_carga)
    new_self = calcular_self(Resname, aob, Num_polares, Num_no_polares)

    delta_elec = new_elec - elec
    delta_self = delta_carga * new_self
        
    # Calcular diferencia de energía
    delta_mc = delta_carga * (ph - Pka_ref) * kb * T * np.log(10) + delta_elec + delta_self
        
    if delta_mc < 0:
        # Si la energía disminuye, aceptar la nueva carga
        for index, (residuo, _) in enumerate(lista_residuos):
            if residuo == residuo_mc:
                lista_residuos[index] = (residuo, new_carga)
                lista_cambios.extend(lista_residuos)  # Agregar elementos directamente
    else:
        # Si la energía aumenta, aceptar con una probabilidad determinada por el criterio de Metropolis
        random_prob = random.uniform(0, 1)
        if random_prob <= np.exp(-delta_mc/(kb*T)):# borre el dividido KbT
            for index, (residuo, _) in enumerate(lista_residuos):
                if residuo == residuo_mc:
                    lista_residuos[index] = (residuo, new_carga)
                    lista_cambios.extend(lista_residuos)  # Agregar elementos directamente
        else:
            # Si no se acepta el cambio, agregar la configuración actual a la lista de cambios
            lista_cambios.extend(lista_residuos)  # Agregar elementos directamente
    
 # Filtrar enlaces_matrix para incluir solo los enlaces que involucren al residuo_mc
    bond_matrix = []
    for _, row in enlaces_matrix.iterrows():
        if row['seq_i'] == residuo_mc or row['seq_j'] == residuo_mc:
            bond_index = int(row['bond_index'])  # Asegurar que bond_index es un entero
            # Obtener las cargas correspondientes de lista_residuos
            carga_i = next(carga for residuo, carga in lista_residuos if residuo == row['seq_i'])
            carga_j = next(carga for residuo, carga in lista_residuos if residuo == row['seq_j'])
            bond_matrix.append([bond_index, carga_i, carga_j])
    
    # Convertir bond_matrix a DataFrame
    bond_matrix_df = pd.DataFrame(bond_matrix, columns=['bond_index', 'carga_i', 'carga_j'])
    #print(data)
    
    # Retornar lista de cambios y la matriz
    return lista_cambios, bond_matrix_df
    
    
    


#############################################
##### MAS FUNCIONES #########################
#############################################

def procesador_de_archivo_con_residuos_cargados(archivo):
    with open(archivo, 'r') as file:
        residues = []
        for line in file:
            parts = line.split()
            residue = int(parts[0])
            charge = float(parts[1])
            residues.append((residue, charge))

    # Filtrar los residuos con carga diferente de cero
    charged_residues = [(residue, charge) for residue, charge in residues if charge != 0.0]
    
    return charged_residues


def constructor_df_de_posiciones_de_CB(cb_atoms_matrix, positions):
    cb_positions = []

    # Iterar sobre los índices de átomos CB y números de residuo en cb_atoms_matrix
    for atom_index, residue_number in cb_atoms_matrix:
        cb_position = positions[atom_index]
        cb_positions.append([residue_number-1, atom_index, cb_position.x, cb_position.y, cb_position.z])
    
    # Definir las columnas del DataFrame
    columns = ['Residue_Number', 'Atom_Index', 'X', 'Y', 'Z']
    cb_positions_df = pd.DataFrame(cb_positions, columns=columns)
    
    return cb_positions_df
    
def constructor_df_de_los_objetos_force(last_force, num_bonds, atom_list, oa):
    enlaces_info = []

    for bond_index in range(num_bonds):
        bond_parameters = last_force.getBondParameters(bond_index)
        atom_index_i = bond_parameters[0]
        atom_index_j = bond_parameters[1]
        charge_i, charge_j = bond_parameters[2]
        residue_index_i = atom_list[atom_index_i].residue.index
        residue_index_j = atom_list[atom_index_j].residue.index
        num_seq_i = oa.residues[residue_index_i].index
        num_seq_j = oa.residues[residue_index_j].index
        enlace_info = [bond_index, num_seq_i, num_seq_j, charge_i, charge_j]
        enlaces_info.append(enlace_info)

    columns = ['bond_index', 'seq_i', 'seq_j', 'carga_i', 'carga_j']
    enlaces_matrix = pd.DataFrame(enlaces_info, columns=columns)
    
    return enlaces_matrix
    
    
def id_and_residues_df_oasistem(oa):
    cb_atoms_matrix = [
        (atom.index, int(residue.id))
        for residue in oa.pdb.topology.residues()
        for atom in residue.atoms()
        if atom.name == 'CB'
    ]
    return cb_atoms_matrix

