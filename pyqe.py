from pandas import DataFrame, read_csv
from numpy import array,arange, sqrt, concatenate
from scipy import interpolate
import os
import subprocess
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import warnings
warnings.filterwarnings('ignore')

####global functions
#--------------------------------------------------------------------------------------------------------------

#general functions

#find function, use to find a word in a file. input reg: string, input linhas: list of list from the orifinal file
def find(reg, linhas):  
    # cria a lista que vai ser retornada como resultado da funcao com a palavra e seu valor
    busca_resultado = []
    # for para varrer todas as linhas, o valor de i indica qual á linha em linhas[i]
    for i in range(len(linhas)):
        # for para varrer cada linha a procura da palavra buscada
        for item in linhas[i]:
            # verifica se o registro esta na linha
            if item == reg:
                # este for e para varrer novamente a linha se achar o registro para gravar todos os valores da linha em uma lista
                for item in linhas[i]:
                    if item != reg:
                        busca_resultado.append(item)
    return busca_resultado

#find n line function, use to find a line number of a word in a file. input reg: string; input linhas: list of list from the orifinal file
def find_n_line(reg, linhas):
    # cria a lista que vai ser retornada como resultado da funcao com a palavra e seu valor
    busca_resultado = []
    # for para varrer todas as linhas, o valor de i indica qual á linha em linhas[i]
    j = 0
    for i in range(0, len(linhas)):
        # for para varrer cada linha a procura da palavra buscada
        for item in linhas[i]:
            # verifica se o registro esta na linha
            if item == reg:
                # este for e para varrer novamente a linha se achar o registro para gravar todos os valores da linha em uma lista
                busca_resultado.append(j)
        j += 1
    return busca_resultado


# function that transform file in to list. input name: path and the name of the file.
def open_file(name):
    file = open(name, 'r', encoding='utf-8',
                errors='ignore')  # open the file.
    # creat a list that will be full filed with the file with a line as a list element.
    file_list = []
    for line in file:  # open a loop that cover the file.
        line = line.strip('\n')  # drop out all '\n' contained in every line.
        # change the spaces for a element of a list, ex: 'the energy is' --> ['the','energy','is'].
        line = line.split()
        file_list.append(line)  # add the line in the list file_list.
    file.close()  # close de file.
    return file_list

def number_identify(string):
        number = ''
        number_list = list(map(lambda x:str(x),range(10)))
        for caracter in string:
            if (caracter in number_list) == True:
                number+=caracter
            elif caracter == '.':
                number+=caracter
            elif caracter == '+':
                old_number = number
                number = ''
        return float(number)

#--------------------------------------------------------------------------------------------------------------

#functions to get information from scf.out, bands.out and nscf.out

#function to get the quantum expresso version from output. input: file, list of list from original output file
def version_quantum(file):
    for line in file:
        for item in line:
            if item == 'Program':
                version = line[2][2:]
    return version

#function to get the fermy energy from outpu. input: file, list of list from original output file; input type: string for whitch ouput, scf.out,nscf.out or bands.out; input version: string with the version from quantum expresso 
def E_fermi(file,type,version):
    i =0
    for line in file:
        for item in line:
            if version == '7.1':
                if type == 'bands.out' or type == 'nscf.out':
                    if item == "Fermi":
                        count = i
                    elif item == "occupied," or item == "occupied":
                        count = i
                elif type == 'scf.out':
                    if item == "Fermi":
                        count = i
                    elif item == "occupied," or item == "occupied":
                        count = i    
            elif version == '7.0':
                if type == 'bands.out' or type == 'nscf.out':
                    if item == "Fermi":
                        count = i
                    elif item == "occupied," or item == "occupied":
                        count = i
                elif type == 'scf.out':
                    if item == "Fermi":
                        count = i
                    elif item == "occupied," or item == "occupied":
                        count = i                      
                        
            elif version == '6.7MaX':
                if type == 'scf.out':
                    if item == "Fermi":
                        count = i
                    elif item == "occupied," or  item ==  "occupied":
                        count = i    
                
        i+=1
    for item in file[count]:
        try:
            E_f = round(float(item), 4)
            break
        except ValueError:
                pass
    return E_f

#function to get the volume of cell. input: file, list of list from original output file
def volumn(file):
    for line in file:
        for i in range(len(line)):
            if line[i] == 'volume':
                vol =  float(line[i+2])
                return vol


#function to get the a latice. input: file, list of list from original output file
def alat(file):
    for line in file:
        for i in range(len(line)):
            if line[0]== "lattice" and line[1] == 'parameter':
                return float(line[4])

#--------------------------------------------------------------------------------------------------------------

#functions for pdos files

#function to correct the head for pdos files. input arquivo , file with the path from the original pdos file
def change_coluns_pdos(arquivo):
    # abre o arquivo scf.out
    arq = open(arquivo, 'r')
    # cria a lista que contenhaa cada linha do arquvio cif como um elemento da lista
    new_arq = []
    ct = 0
    for linha in arq:
        # retira o \n de quebra de linnha da linha
        line = linha.strip('\n')
        # retira os espaços entre os elementos da linha, assim criando uma lista
        line = line.split()
        if ct == 0 and line[0] == '#':
            line.remove('#')
            new_line = []
            for item in line:
                if item != '(eV)':
                    new_line.append(item)
        else:
            new_line = line
        new_arq.append(new_line)
        ct += 1
    # fecha o arquivo
    arq.close()

    arq = open(arquivo,'w')
    for item in new_arq:
        line = ''
        for word in item:
            line+= word + '   '
        line += '\n'
        arq.writelines(line)
    arq.close()
    return 0

#function to get the atoms from scf.in file; input scf_in_file: list os list from original scf.in file
def atoms(scf_in_file):
    i = 0;atom = []
    for line in scf_in_file:
        for item in line:
            if item == 'ATOMIC_SPECIES' or item == 'atomic_species':
                count1 = i+1
            elif item == 'ATOMIC_POSITIONS' or item == 'atomic_positions':
                count2 = i
        i+=1
    for j in range(count1,count2):
        atom.append(scf_in_file[j][0])
    return atom


#function to get all of pdos files and covert in to a dicitionary delimted for orbiatis( s,p,etc) and for atoms. input: path , path from the pdos files are; input prefix, the prefix if exist in the pdos files; input atoms, atoms get from the functions atoms
def pdos_files_get(path,prefix,atoms):
    pdos_files = []
    for dir,subdir,files in os.walk(path):
        for file in files:
            if file[:len(prefix)] == prefix and file != prefix+'.pdos.dat.pdos_tot' and file != prefix + '.pdos-proj.dat.projwfc_up' and file != prefix+'.pdos_tot':
                pdos_files.append(os.path.join(path,file))

    pdos_atom_dict = {};pdos_orb_dict = {}
    for a in range(len(atoms)):
        list = [];orb_list = []
        for name in pdos_files:
            for i in range(1,len(name)):
                if len(atoms[a]) == 1:
                    if name[i] == '(' and name[i+1] == atoms[a]:
                        list.append(name)
                    elif name[i] == '(' and name[i+1] == 's':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'p':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'd':
                       if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'f':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])                       
                else:
                    if name[i-1] == '(' and name[i] == atoms[a][0] and name[i+1] == atoms[a][1]:
                        list.append(name)
                    elif name[i] == '(' and name[i+1] == 's':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'p':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'd':
                       if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])
                    elif name[i] == '(' and name[i+1] == 'f':
                        if (name[i+1] in orb_list) == False:
                            orb_list.append(name[i+1])     

        pdos_atom_dict[atoms[a]] = list
        pdos_orb_dict[atoms[a]] = orb_list

    pdos_atom_orb_dict = {}
    for atom in atoms:
        pdos_atom_orb_dict[atom] = {}
        for orb in pdos_orb_dict[atom]:
            pdos_atom_orb_dict[atom][orb] = []
            list = []
            for file in pdos_atom_dict[atom]:
                for i in range(len(file)):
                    if file[i] == '(' and file[i+1] == orb:
                        list.append(file)
            pdos_atom_orb_dict[atom][orb] = list

    return [pdos_atom_dict,pdos_atom_orb_dict]


#--------------------------------------------------------------------------------------------------------------

#functions for bands files
#function to transform the bands.dat.gnu file in to a data frame. input arquivo, is the orifinal aquivo bands.dat.gnu
def bandas_df(arquivo):
    bandas = read_csv(arquivo, delim_whitespace=True,
                         names=['k', 'E'], dtype=float)
    # subtrai todos os valores pelo nivel de fermi
    return bandas

#function to transform the original dataframe from function bandas_df for a new dataframe separeted for orbitals. input bandas, dataframme obtained from function bandas_df; input ef, bool varable, if u want to correct the fermi level or not; input Ef: float, the value of fermyenergy
def orbitais(bandas, ef, Ef):
    # nesse bloco sera feito a separaçao do data frame em varias colunas, onde cada coluna é referente a uma banda
    # neste passo vamos obter do dataFrame Bandas original apenas os valores de momento que nao são repetidos
    # cria a lista para armazenaros valroes
    if ef == True:
        bandas['E'] -= Ef
    momentos = []
    # for para varrer todos os indices do dataframe
    for i in bandas.index:
        # pega o valor do momento para o indice i, no caso a linha i do data frame
        mom = bandas.loc[i][0]
        # cria uma condicao para adicionar na lista apenas itens nao repetidos, ou seja se mom nao estiver na  lista momentos, entao ele sera adicionado a mesma, caso contrario nada e feito
        condicao = mom in momentos
        if condicao == False:
            momentos.append(mom)
    # variavel para fazer a contagem total do numero de bandas
    n_bnd = 0
    # para fazer a contagem vamso pegar todos os valores de k na coluna "k" do dataframe Bandas, sempre que k=0, quer dize que começamos uma nova banda, portanto a
    # condicao de troca de bandas sera quando k=0, sempre que isso acontecer, adicionamos +1 no valor de n_bnd(numero de bandas)
    cond = bandas["k"] == 0
    # começa o for varrendo todos os valores de k contidos na coluna bandas["k"]
    for item in cond:
        if item == True:
            n_bnd += 1
        else:
            pass
    # cria o novo dataFrame com a coluna contendo apenas os valores de k que nao sao repetidos
    df = DataFrame(momentos, columns=["K"], dtype=float)
    # começa com o menor valor posivel para a energia da banda de valencia
    E_valencia = min(bandas['E'])
    # começa com o maior valor posivel para a energia da banda de conduçao
    E_conducao = max(bandas['E'])
    gap = True
    # loop para varrer todas as bandas do data frame
    for j in range(n_bnd):
        # cria a lista que vai armazear os valores de energia para a banda j
        energia = []
        # loop para varrer todos os valores de energia contidos no dataframe original que coresponem a banda j
        for i in range(len(momentos)):
            e = bandas.loc[i+j*len(momentos)][1]
            energia.append(e)
        # cria uma nova coluna no novo dataframe adicionando todos os valores de energia da banda j
        df["orb "+str(j)] = energia
        #verifica se é metal
        e_min = min(energia); e_max = max(energia)
        if ef == False:
            if e_min <Ef and e_max>Ef:
                gap = False
        else:
            if e_min <0 and e_max>0:
                gap = False

        # atualiza o valor da energia de valencia
        if ef == False:
            if max(energia) >= E_valencia and max(energia) <= Ef:
                E_valencia = max(energia)
                for k in range(len(energia)):
                    if energia[k] == E_valencia:
                        E_valencia_point = [momentos[k],energia[k]]
            # atualiza o valor da energia de conduçao
            elif min(energia) <= E_conducao and min(energia) >= Ef:
                E_conducao = min(energia)
                for k in range(len(energia)):
                    if energia[k] == E_conducao:
                        E_conducao_point = [momentos[k],energia[k]]
        else:
            if max(energia) >= E_valencia and max(energia) <= 0:
                E_valencia = max(energia)
                for k in range(len(energia)):
                    if energia[k] == E_valencia:
                        E_valencia_point = [momentos[k],energia[k]]                
            # atualiza o valor da energia de conduçao
            elif min(energia) <= E_conducao and min(energia) >= 0:
                E_conducao = min(energia)
                for k in range(len(energia)):
                    if energia[k] == E_conducao:
                        E_conducao_point = [momentos[k],energia[k]]
    if gap == True:
        return [df, round(E_valencia,3), round(E_conducao,3),E_valencia_point,E_conducao_point]
    else:
        if ef == False:
            return [df, Ef, Ef,0,0]
        else:
            return [df,0, 0,0,0]

#function to get the kpoints letter from the bands.in. input bandin_file, list of list from the original bands.in file
def k_point(bandin_file):
    k_points_band = []
    k_points_letter = []
    nk = []
    i = 0
    for line in bandin_file:
        if len(line) >= 1:
            if line[0] == 'K_POINTS':
                n_line = i
        i += 1
    n_kpoints = int(bandin_file[n_line+1][0])
    for i in range(n_line+2, n_line+2+n_kpoints):
        k_points_band.append(array(bandin_file[i][:3], dtype=float))
        nk.append(int(bandin_file[i][3]))
        if bandin_file[i][4][1:] == 'G':
            k_points_letter.append('$\Gamma$')
        else:
            k_points_letter.append(bandin_file[i][4][1:])
    return [k_points_letter, k_points_band,nk]

#function to que n_spin from bands.in.
def n_spin(bandsin_file):
    marker = False
    for line in bandsin_file:
        if len(line) >=1:
            if line[0] == 'nspin=':
                marker = True
                nspin = int(line[1])
            elif line[0] == 'nspin':
                marker = True
                nspin = int(line[2])
            elif line[0][:-1] == 'nspin=':
                marker = True
                nspin = int(line[0][-1:])
                
    if marker == False:
        return 1
    else:
        return nspin

def verbosity(bandsin_file):
    marker = False
    for line in bandsin_file:
        if len(line) >=1:
            if line[0] == 'verbosity=':
                marker = True
                verbo = line[1][1:-1]
            elif line[0] == 'verbosity':
                marker = True
                verbo = line[2][1:-1]
            elif line[0][:-1] == 'verbosity=':
                marker = True
                verbo= line[0][-1:][1:-1]
    if marker == False:
        return 'low'
    else:
        return verbo
#funtion to get the momentum axis values from bands.x.out file. input file, list of list from the original bands.x.out file. input n_k, number of high simmetric points in brilluion zone, geted from bands.in file.
def k_points_path(file, n_k):
    i = 0
    k_points_bandsx = []
    k_path = []
    for line in file:
        for item in line:
            if item == 'wavefunctions':
                n_line = i+1
        i += 1
    for j in range(n_line, n_line + n_k):
        
        
        try:
            n1 = file[j][2]
            n2 = file[j][3]
            n3 = file[j][4]
            k_points_bandsx.append(array(file[j][2:5], dtype=float))
            k_path.append(float(file[j][7]))
        except ValueError:
            try:
                n1 = float( file[j][2] )
                for k in range(len(file[j][3])):
                    if file[j][3][k] == '-':
                        count = k 
                        break
                n2 = float(file[j][3][:count])
                n3 = float(file[j][3][count:])
            except ValueError:
                for k in range(1,len(file[j][2])):
                    if file[j][2][k] == '-':
                        count = k 
                        break
                n1 = float(file[j][2][:count])
                try:
                    n2 = float(file[j][2][count:])
                    n3  = float(file[j][3])
                    
                except ValueError:
                    for k in range(count,len(file[j][2])):
                        if file[j][2][k] == '-':
                            count2 = k 
                            break
                        n2 = float(file[j][2][count:count2])
                        n3  = float(file[j][2][count2:])
            k_points_bandsx.append([n1,n2,n3] )
            l = 0
            for item in file[j]:
                if item == 'coordinate':
                    k_path.append(float(file[j][l+1]))
                    break
                l+=1
    return [k_path, k_points_bandsx]


#--------------------------------------------------------------------------------------------------------------
#functions for projected bands

#function to get projeceted contribuitions from every state on atomitc_proj.xml. Input: File, list of list from the original file
def atomic_proj(file):
    def number_identify(string):
        number = ''
        number_list = list(map(lambda x:str(x),range(10)))
        for caracter in string:
            if (caracter in number_list) == True:
                number+=caracter
            elif caracter == '.':
                number+=caracter
            elif caracter == '+':
                old_number = number
                number = ''
        return float(number)

    for line in file:
        if line[0] == '<HEADER':
            n_band = int(number_identify(line[1]))
            n_spin = int(number_identify(line[3]))
            n_wfc = int(number_identify(line[4]))


    n = 0
    for i in range(len(file)):
            if file[i][0] == '<PROJS>':
                ct = i
            elif file[i][0] == '</PROJS>':
                    n = i-ct
            elif n != 0:
                break
    wfc_proj = {}
    for i in range(1,n_wfc+1):
        wfc_proj[i] = {}
        for j in range(1,n_band+1):
            wfc_proj[i][j] = []
    for wfc in range(1,n_wfc+1):
            for i in range(len(file)):
                if file[i][0] == '<ATOMIC_WFC':
                    #print(file[i])
                    index = int(number_identify(file[i][1]))
                    if index == wfc:
                        ct = 1
                        for j in range(i+1,i+1+n_band):
                            psi_real = float(file[j][0]);psi_imag = float(file[j][1])
                            psi_2 = pow(psi_real,2) +pow(psi_imag,2)
                            wfc_proj[wfc][ct].append(psi_2)
                            ct+=1          
    return [wfc_proj,n_band,n_spin,n_wfc]

#function to get states from bn.prowfc.out. Input:File, list of list from the original file
def States(file):
    states_line = []; atom_number = {};atom_number2 = {}
    for line in file:
        if len(line) >=1:
            if line[0] == 'state':
                states_line.append(line)

    for line in states_line:
        atom_number2[line[5][1:]] = []
    for line in states_line:
        if (int(line[4])in atom_number2[line[5][1:]]) == False:
            atom_number2[line[5][1:]].append(int(line[4]))
        atom_number[int(line[4])] = line[5][1:]
    atom_orb = {}
    for key in atom_number.keys():
        marc_s = False; marc_p = False; marc_d = False
        for line in states_line:
            if int(line[4]) == key:
                if line[8] == '1':
                    if marc_s == False:
                        atom_orb[key] = {'s':[line[2][:-1]]}
                        marc_s = True
                    else:
                        atom_orb[key]['s'].append(line[2][:-1])
                elif line[8] == '2':
                    if marc_p == False:
                        atom_orb[key].update({'p':[line[2][:-1]]})
                        marc_p = True
                    else:
                        atom_orb[key]['p'].append(line[2][:-1])
                elif line[8] == '3':
                    if marc_d == False:
                        atom_orb[key].update({'d':[line[2][:-1]]})
                        marc_d = True
                    else:
                        atom_orb[key]['d'].append(line[2][:-1])
        
    return [atom_number2,atom_orb]

#function to transform all information in to dictonary. Input: atom_number, dict with atom and the correspondent wfc number;Input: atom_number_orb, dict with atom and the correspondent state, Input: proj, geted from atomic_proj function, Input: nbands, number of bands
def projcted_atoms(atom_number,atom_orb_number,proj,n_band):
    atom_orb = {}
    for atom in atom_number.keys():
        atom_orb[atom] = {}
        for number in atom_number[atom]:
            for orb in atom_orb_number[number].keys():
                s=0;p=0;d=0;marker_s = False;marker_p = False;marker_d = False
                for item in atom_orb_number[number][orb]:
                    if orb == 's':                      
                        atom_orb[atom][orb] = proj[int(item)]
                    elif orb == 'p':
                        if marker_p == False:
                            atom_orb[atom][orb] = {}
                            marker_p = True
                        if p == 0:
                            atom_orb[atom][orb].update({'pz':proj[int(item)]})
                        elif p == 1:
                            atom_orb[atom][orb].update({'px':proj[int(item)]})
                        elif p == 2:
                            atom_orb[atom][orb].update({'py':proj[int(item)]})
                        p+=1
                    elif orb == 'd':
                        if marker_d == False:
                            atom_orb[atom][orb] = {}
                            marker_d = True
                        if d == 0:
                            atom_orb[atom][orb].update({'dz2':proj[int(item)]})
                        elif d == 1:
                            atom_orb[atom][orb].update({'dzx':proj[int(item)]})
                        elif d == 2:
                            atom_orb[atom][orb].update({'dzy':proj[int(item)]})
                        elif d == 3:
                            atom_orb[atom][orb].update({'dx2-y2':proj[int(item)]})
                        elif d == 4:
                            atom_orb[atom][orb].update({'dxy':proj[int(item)]})
                        d+=1                        

    for atom in atom_orb.keys():    
        temp = {};temp2 = {}; temp3 = {}
        for i in range(1,n_band+1):
            for orb in atom_orb[atom]:
                marker = False
                if orb == 's':
                    temp2[i] = array(atom_orb[atom][orb][i])
                    temp3[i] = array(atom_orb[atom][orb][i])
                elif orb == 'p':
                    for item in atom_orb[atom][orb]:
                        if marker == False:
                            temp[i] = array(atom_orb[atom][orb][item][i])
                            marker = True
                        else:
                            temp[i] = temp[i] + array(atom_orb[atom][orb][item][i])
                        temp2[i] = temp2[i] + temp[i]
                        temp3[i] = temp3[i] + temp[i]
                elif orb == 'd':
                    for item in atom_orb[atom][orb]:
                        if marker == False:
                            temp[i] = array(atom_orb[atom][orb][item][i])
                            marker = True
                        else:
                            temp[i] = temp[i] + array(atom_orb[atom][orb][item][i])
                        temp2[i] = temp2[i] + temp[i]
        atom_orb[atom].update({'s+p':temp3})                        
        atom_orb[atom].update({'tot':temp2})            
        atom_orb[atom]['p'].update({'tot':temp})
    return atom_orb

#function to get band data from bands.out. Input: file, list of list from the original bands.out file
def band_structure_from_bandsout(file, fermy_energ_cor, fermy_energy,version):
    n_band = 0
    for line in file:
        if len(line)>=3:
            if line[2] == 'Kohn-Sham':
                n_band = int(line[4])
    bands = {}
    for e in range(1,n_band+1):
        bands[e] = []
    k_value = [];k_tol = 0
    for i in range(len(file)):
        if len(file[i])>=1:
            if file[i][0] == 'k':
                line = file[i]
                if len(k_value)>=1:
                    if ( '-' in line[1]) == True:
                        try:
                            kx = float(line[1][1:])
                            ky = float(line[2])
                            kz = float(line[3])
                        except ValueError:
                            try:
                                kx = float(line[1][1:])
                                for k in range(len(line[2])):
                                    if line[2][k] == '-':
                                        count = k 
                                        break
                                ky = float(line[2][:count])
                                kz = float(line[2][count:])
                            except ValueError:
                                for k in range(2,len(line[1])):
                                    if line[1][k] == '-':
                                        count = k 
                                        break
                                kx = float(line[1][1:count])
                                try:
                                    ky = float(line[1][count:])
                                    kz  = float(line[2])
                                except ValueError:
                                    for k in range(count+1,len(line[1])):
                                        if line[1][k] == '-':
                                            count2 = k 
                                            break
                                    ky = float(line[1][count:count2])
                                    kz  = float(line[1][count2:])
                    
                    else:
                    
                        try:
                            kx = float(line[2])
                            ky =  float(line[3])
                            kz = float(line[4])
                        except ValueError:
                            try:
                                kx = float(line[2])
                                for k in range(len(line[3])):
                                    if line[3][k] == '-':
                                        count = k 
                                        break
                                ky = float(line[3][:count])
                                kz = float(line[3][count:])
                            except ValueError:
                                for k in range(1,len(line[2])):
                                    if line[2][k] == '-':
                                        count = k 
                                        break
                                kx = float(line[2][:count])
                                try:
                                    ky = float(line[2][count:])
                                    kz  = float(line[3])
                                    
                                except ValueError:
                                    for k in range(count+1,len(line[2])):
                                        if line[2][k] == '-':
                                            count2 = k 
                                            break
                                    ky = float(line[2][count:count2])
                                    kz  = float(line[2][count2:])
                         
    
                    k_mod = sqrt((kx-kx_old)**2+(ky-ky_old)**2+(kz-kz_old**2))
                    k_tol+=k_mod

                    k_value.append(round(k_tol,4))
                    del(kx_old);del(ky_old);del(kz_old)
                    kx_old = kx; ky_old = ky; kz_old = kz
                else:
                    try:
                        kx = float(line[2])
                        ky =  float(line[3])
                        kz = float(line[4])
                    except ValueError:
                        try:
                            kx = float(line[2])
                            for k in range(len(line[3])):
                                if line[3][k] == '-':
                                    count = k 
                                    break
                            ky = float(line[3][:count])
                            kz = float(line[3][count:])
                        except ValueError:
                            for k in range(1,len(line[2])):
                                if line[2][k] == '-':
                                    count = k 
                                    break
                            kx = float(line[2][:count])
                            try:
                                ky = float(line[2][count:])
                                kz  = float(line[3])
                                
                            except ValueError:
                                for k in range(count,len(line[2])):
                                    if line[2][k] == '-':
                                        count2 = k 
                                        break
                                    ky = float(line[2][count:count2])
                                    kz  = float(line[2][count2:])
                         
            

                    k_mod = sqrt(kx**2+ky**2+kz**2)
                    k_tol+=k_mod
                    k_value.append(round(k_tol,4))
                    kx_old = kx; ky_old = ky; kz_old = kz
                
                bands_line = []
                for k in range(i+1,i+len(file)):
                    if len(file[k]) >=1: 
                        if file[k][0] == 'k' or file[k][0] == 'Writing' or file[k][0] == 'highest':
                            break
                        else:
                            if version == '7.1' or version == '7.0':
                                if file[k][0] != ' ':
                                    for item in file[k]:
                                        bands_line.append(float(item))
                            elif version == '6.7MaX' or version == '6.7':
                                if file[k][0] != ' ' and file[k][0] != 'occupation' and file[k][0] != '1.0000' and file[k][0] != '0.0000':
                                    for item in file[k]:
                                        bands_line.append(float(item))
                for j in range(n_band):
                    bands[j+1].append(bands_line[j])
    df = DataFrame({'k':k_value})
    for e in range(1,n_band+1):
        if fermy_energ_cor == True:
            df[str(e)]= array(bands[e]) - fermy_energy
        elif fermy_energ_cor == False:
            df[str(e)]= bands[e] 
        else:
            print('Worng value for fermy energy correction parameter, only True and False are accepted. Default: True')
    return df

#--------------------------------------------------------------------------------------------------------------
#functions for change density

def get_vectors(file):
    e11 = 'e1(1)';e12 = 'e1(2)';e13 = 'e1(3)'
    e1 = [e11,e12,e13]
    e21 = 'e2(1)';e22 = 'e2(2)';e23 = 'e2(3)'
    e2 = [e21,e22,e23]
    x1 = 'x0(1)';x2 = 'x0(2)';x3= 'x0(3)'
    x_list = [x1,x2,x3]
    e1_v = []; e2_v = []; x_0 = []
    for e in e1:
        for line in file:
            for i in range(len(line)):
                if line[i] == e:
                    if (',' in line[i+2]):
                        e1_v.append(float(line[i+2][:-1]))
                    else:                  
                        e1_v.append(float(line[i+2]))
                        
    for e in e2:
        for line in file:
            for i in range(len(line)):
                if line[i] == e:
                    if (',' in line[i+2]):
                        e2_v.append(float(line[i+2][:-1]))
                    else:                  
                        e2_v.append(float(line[i+2]))
                        
    for x in x_list:
        for line in file:
            for i in range(len(line)):
                if line[i] == x:
                    if (',' in line[i+2]):
                        x_0.append(float(line[i+2][:-1]))
                    else:                  
                        x_0.append(float(line[i+2]))
                        
                        
    for line in file:
       for i in range(len(line)):
           if line[i] == 'nx':
               nx = line[i+2]
           elif line[i] == 'ny':
               ny = line[i+2]
    return [e1_v,e2_v,x_0,int(nx),int(ny)]


#classes
class File:  # this class is relative to treatement of the raw data in to a better analize format
    def __init__(self):  # load the input file
        self.file_dir = './'
        self.scfin_file_name = 'scf.in'
        self.scfout_file_name = 'scf.out'
        # must have ! in front of k letter, like 0.3333333 0.3333333 0.00000000 !K , for gamma point use G
        self.bandin_file_name = 'bands.in'
        self.bandout_file_name = 'bands.out'
        self.bandxout_file_name = 'bands.x.out'
        self.band_file_name = 'bands.dat.gnu'
        self.dos_file_name = 'dos'
        self.pdos_file_name = 'pdos.dat.pdos_tot'
        self.pdos_prefix = ''
        self.fermi_energy = 0
        self.projwfcout_file_name = 'projwfc.out'
        self.atomic_xml_name = 'atomic_proj.xml'
        self.ppin_file_name = 'pp.in'
        self.charge2D = 'charge2D.out'

    def reset_atributes(self):
        del(self.file_dir)
        del(self.scfin_file_name)
        del(self.scfout_file_name)
        del(self.bandin_file_name)
        del(self.bandout_file_name)
        del(self.bandxout_file_name)
        del(self.band_file_name)
        del(self.dos_file_name)
        del(self.pdos_file_name)
        del(self.pdos_prefix)
        del(self.fermi_energy)
        del(self.projwfcout_file_name )
        del(self.atomic_xml_name)
        del(self.ppin_file_name)
        del(self.charge2D)


        self.file_dir = './'
        self.scfin_file_name = 'scf.in'
        self.scfout_file_name = 'scf.out'
        self.nscfin_file_name = 'nscf.in'
        self.nscfout_file_name = 'nscf.out'
        # must have ! in front of k letter, like 0.3333333 0.3333333 0.00000000 !K , for gamma point use G
        self.bandin_file_name = 'bands.in'
        self.bandout_file_name = 'bands.out'
        self.bandxout_file_name = 'bands.x.out'
        self.band_file_name = 'bands.dat.gnu'
        self.dos_file_name = 'dos'
        self.pdos_file_name = 'pdos.dat.pdos_tot'
        self.pdos_prefix = ''
        self.fermi_energy = 0
        self.projwfcout_file_name = 'projwfc.out'
        self.atomic_xml_name = 'atomic_proj.xml'
        self.ppin_file_name = 'pp.in'
        self.charge2D = 'charge2D.out'      
        


    def files_atributes(self):
        print("Files path: files_folder; scf.in: scf_in_name; scf.out: scf_out_name; nscf.in: nscf_in_name; nscf.out: nscf_out_name; bands.in: bands_in_name; bands.out: bands_out_name; bands.x.out: bandsx_name; band_dat.gnu: bands_dat_name, dos.dat: dos_dat_name; pdos.dat: pdos_dat_name; pdos prefix: pdos_prefix; bandsup.dat.gnu: bandsup_dat_name; bandsdown.dat.gnu: bandsdown_dat_name; projwfc.out: projwfc_out_name; ppin_name: pp.in file name; charge2D_dat_name: charge2D.out "  )

    def Set_files_atributes(self, dict):
        keys = dict.keys()
        for key in keys:
            if key == 'files_folder':
                self.file_dir = dict[key]
            elif key == 'scf_in_name':
                self.scfin_file_name = dict[key]
            elif key == 'scf_out_name':
                self.scfout_file_name = dict[key]
            elif key == 'nscf_in_name':
                self.nscfin_file_name = dict[key]
            elif key == 'nscf_out_name':
                self.nscfout_file_name = dict[key]
            elif key == 'bands_in_name':
                self.bandin_file_name = dict[key]
            elif key == 'bands_out_name':
                self.bandout_file_name = dict[key]
            elif key == 'bandsx_name':
                self.bandxout_file_name = dict[key]
            elif key == 'bands_dat_name':
                self.band_file_name = dict[key]
            elif key == 'dos_dat_name':
                self.dos_file_name = dict[key]
            elif key == 'pdos_dat_name':
                self.pdos_file_name = dict[key]
            elif key == 'pdos_prefix':
                self.pdos_prefix = dict[key]
            elif key == 'bandsup_dat_name':
                self.bandup_file_name = dict[key]
            elif key == 'bandsdown_dat_name':
                self.banddown_file_name = dict[key]
            elif key == 'projwfc_out_name':
                self.projwfcout_file_name = dict[key]
            elif key == 'ppin_name':
                self.ppin_file_name = dict[key]
            elif key == 'charge2D_dat_name':
                self.charge2D_name = dict[key]
                

    def Load(self):
        self.scfin_file_path = os.path.join(
            self.file_dir, self.scfin_file_name)
        self.scfin_file = open_file(self.scfin_file_path)

        self.scfout_file_path = os.path.join(
            self.file_dir, self.scfout_file_name)
        self.scfout_file = open_file(self.scfout_file_path)
        self.version_QE = version_quantum(self.scfout_file)
        self.fermi_energy = E_fermi(self.scfout_file,'scf.out',self.version_QE)
        self.cell_volumn = volumn(self.scfout_file)
        self.alat = alat(self.scfout_file)

    def Bands_files(self, fermi_energy_corr=True, bands_project = False):
        self.bandout_file_path = os.path.join(self.file_dir, self.bandout_file_name)
        self.bandout_file = open_file(self.bandout_file_path)

        self.bandin_file_path = os.path.join(
            self.file_dir, self.bandin_file_name)
        self.bandin_file = open_file(self.bandin_file_path)
        self.nspin = n_spin(self.bandin_file)
        self.k_points_letter, self.k_points_bands, self.nk = k_point(self.bandin_file)

        self.verbosity = verbosity(self.bandin_file)
            
        self.version_QE = version_quantum(self.bandout_file)
        
        if self.version_QE == '7.1' and self.verbosity == 'low':
            self.fermi_energy = E_fermi(self.bandout_file,'bands.out',self.version_QE)


        if self.nspin == 1:
            self.band_file_path = os.path.join(self.file_dir, self.band_file_name)
            bands = bandas_df(self.band_file_path)
            self.bands, self.E_valence, self.E_conduction, self.E_valence_point,self.E_conduction_point = orbitais(
            bands, fermi_energy_corr, self.fermi_energy)
            self.gap = self.E_conduction - self.E_valence

        elif self.nspin == 2:
            self.bandup_file_path = os.path.join(self.file_dir, self.bandup_file_name)
            bands_up = bandas_df(self.bandup_file_path)
            self.bands_up, self.E_valence_up, self.E_conduction_up, self.E_valence_point_up,self.E_conduction_point_up = orbitais(
                bands_up, fermi_energy_corr, self.fermi_energy)

            self.banddown_file_path = os.path.join(self.file_dir, self.banddown_file_name)
            bands_down = bandas_df(self.banddown_file_path)
            self.bands_down, self.E_valence_down, self.E_conduction_down, self.E_valence_point_down,self.E_conduction_point_down = orbitais(
                bands_down, fermi_energy_corr, self.fermi_energy)

            columns = self.bands_up.columns
            self.bands = DataFrame(self.bands_up[columns[0]])
            for i in range(1,len(columns)):
                self.bands[columns[i]] = self.bands_up[columns[i]] + self.bands_down[columns[i]]

            if self.E_conduction_down <= self.E_conduction_up:
                self.E_conduction = self.E_conduction_down
                self.E_conduction_point = self.E_conduction_down
            else:
                self.E_conduction = self.E_conduction_up
                self.E_conduction_point = self.E_conduction_up
            
            if self.E_valence_down >= self.E_valence_up:
                self.E_valence= self.E_valence_down
                self.E_valence_point = self.E_valence_down
            else:
                self.E_valence = self.E_valence_up
                self.E_valence_point = self.E_valence_up


            self.gap = self.E_conduction - self.E_valence
            
        else:
            print('Worng value for nspin atribute.')

        self.bandxout_file_path = os.path.join(
                self.file_dir, self.bandxout_file_name)
        self.bandxout_file = open_file(self.bandxout_file_path)
        self.k_path, self.k_points_bandsx = k_points_path(
                self.bandxout_file, len(self.k_points_letter))

    def Dos_files(self, calculation, fermi_energy_corr=True):  # dos via projwfc.x, nao usar dos.x
        self.nscfout_file_path = os.path.join(self.file_dir, self.nscfout_file_name)
        self.nscfout_file = open_file(self.nscfout_file_path)
        self.version_QE = version_quantum(self.nscfout_file)
        if self.version_QE == '7.1':
            self.fermi_energy = E_fermi(self.nscfout_file,'nscf.out',self.version_QE)

        if calculation == 'projwfc.x':
            self.dos_file_path = os.path.join(
                self.file_dir, self.pdos_file_name)
            change_coluns_pdos(self.dos_file_path)
            self.dos_file = open_file(self.dos_file_path)
            self.dos = DataFrame(self.dos_file)
            self.dos = self.dos.drop(0)
            self.dos = self.dos.rename(columns={0: 'E', 1: 'dos', 2: 'pdos'})
            self.dos = self.dos.astype(float)
            if fermi_energy_corr == True:
                self.dos[self.dos.columns[0]] -= self.fermi_energy
            elif fermi_energy_corr == False:
                pass
            else:
                print(
                    'Wrong value for fermi_energy_corr parameter, only True or False values are accepted. Default = True')

        elif calculation == 'dos.x':
            self.dos_file_path = os.path.join(
                self.file_dir, self.dos_file_name)
            change_coluns_pdos(self.dos_file_path)
            self.dos_file = open_file(self.dos_file_path)
            self.dos_file.remove(self.dos_file[0])
            self.dos = DataFrame(self.dos_file)
            self.dos = self.dos.rename(
                columns={0: 'E', 1: 'dos', 2: 'int_dos'})
            self.dos = self.dos.astype(float)
            if fermi_energy_corr == True:
                self.dos[self.dos.columns[0]] -= self.fermi_energy
            elif fermi_energy_corr == False:
                pass
            else:
                print(
                    'Wrong value for fermi_energy_corr parameter, only True or False values are accepted. Default = True')
        else:
            print(
                "Wrong value for calculation parameter, only 'dos.x' or 'projwfc.x' are accepted.")

    def Pdos_files(self, fermi_energy_corr=True): 

        self.pdos_per_atoms_files, self.pdos_per_orb_files = pdos_files_get(self.file_dir,self.pdos_prefix,atoms(self.scfin_file))  

        self.nscfout_file_path = os.path.join(self.file_dir, self.nscfout_file_name)
        self.nscfout_file = open_file(self.nscfout_file_path)
        self.version_QE = version_quantum(self.nscfout_file)
        if self.version_QE == '7.1':
            self.fermi_energy = E_fermi(self.nscfout_file,'nscf.out',self.version_QE)
        
        #per atoms 
        self.projwfc_out_file_path = os.path.join(self.file_dir, self.projwfcout_file_name)      
        self.projwfc_out_file = open_file(self.projwfc_out_file_path)
        self.atom_number , self.atom_orb_number = States(self.projwfc_out_file)

        pdos_files = {}
        for atom in self.pdos_per_atoms_files.keys():
            i = 0
            for file in self.pdos_per_atoms_files[atom]:
                change_coluns_pdos(file)
                if i ==0:
                    df = read_csv(file, delim_whitespace= True,dtype=float)
                    df.drop(df.columns[2:], axis='columns', inplace=True)
                    maxi = df[df.columns[1]].max()
                else:
                    df2 = read_csv(file, delim_whitespace= True,dtype=float )
                    df[df.columns[1]]+=df2[df2.columns[1]]
                    del(df2)
                i+=1
            if fermi_energy_corr == True:
                df[df.columns[0]] -= self.fermi_energy
                pdos_files[atom] = df
            elif fermi_energy_corr == False:
                pdos_files[atom] = df
            else:
                print('Wrong value for fermi_energy_corr parameter, only True or False values are accepted. Default = True')
        
        self.pdos_per_atoms = pdos_files
        
        #per orb
        orb_files = {}
        for atom in self.pdos_per_orb_files.keys():
            orb_files[atom] = {}
            for orb in self.pdos_per_orb_files[atom]:
                i =0
                for file in self.pdos_per_orb_files[atom][orb]:
                    change_coluns_pdos(file)
                    if i ==0:
                        df = read_csv(file, delim_whitespace= True,dtype=float)
                        df.drop(df.columns[2:], axis='columns', inplace=True)
                    else:
                        df2 = read_csv(file, delim_whitespace= True,dtype=float )
                        df[df.columns[1]]+=df2[df2.columns[1]]
                        del(df2)
                    i+=1
                if fermi_energy_corr == True:
                    df[df.columns[0]] -= self.fermi_energy
                    orb_files[atom][orb] = df

                elif fermi_energy_corr == False:
                    orb_files[atom][orb] = df
                else:
                    print(
                            'Wrong value for fermi_energy_corr parameter, only True or False values are accepted. Default = True')   
        self.pdos_per_orb = orb_files 

    def Bands_pdos(self,fermi_energy_corr = True):
        self.bandout_file_path = os.path.join(self.file_dir, self.bandout_file_name)
        self.bandout_file = open_file(self.bandout_file_path)

        self.bandin_file_path = os.path.join(
            self.file_dir, self.bandin_file_name)
        self.bandin_file = open_file(self.bandin_file_path)
        self.nspin = n_spin(self.bandin_file)
        self.k_points_letter, self.k_points_bands, self.nk = k_point(self.bandin_file)

        self.verbosity = verbosity(self.bandin_file)
            
        self.version_QE = version_quantum(self.bandout_file)
        
        if self.version_QE == '7.1' and self.verbosity == 'low':
            self.fermi_energy = E_fermi(self.bandout_file,'bands.out',self.version_QE)

        self.atomic_xml_file_path = os.path.join(self.file_dir, self.atomic_xml_name)
        self.atomic_xml_file = open_file(self.atomic_xml_file_path)
        self.proj,n_band,nspin,n_wfc = atomic_proj(self.atomic_xml_file)
        

        self.projwfc_out_file_path = os.path.join(self.file_dir, self.projwfcout_file_name)      
        self.projwfc_out_file = open_file(self.projwfc_out_file_path)
        self.atom_number , self.atom_orb_number = States(self.projwfc_out_file)


        #self.band_projec_dict = projcted_atoms(self.atom_number,self.atom_orb_number,self.proj,n_band)
        

        
        if self.verbosity == 'high':
            self.band_df= band_structure_from_bandsout(self.bandout_file, fermi_energy_corr, self.fermi_energy , self.version_QE)

            points_break = [0]
            ct = 0
            for n in self.nk:
                ct +=n
                points_break.append(ct)
            k_axis = []
            k = self.band_df[self.band_df.columns[0]]
            for i in range(len(self.band_df[self.band_df.columns[0]])):
                if (i in points_break) == True:
                    k_axis.append(k[i] )
            self.k_points_path = array(k_axis)
            self.bands_proj = [self.band_df,self.proj,self.k_points_letter,self.k_points_path]
        else:
            print('verbosity must be high for projected bands.')
            
    def Pdos_kresolved(self,fermi_energy_corr = True):
        self.projwfc_out_file_path = os.path.join(self.file_dir, self.projwfcout_file_name)      
        self.projwfc_out_file = open_file(self.projwfc_out_file_path)
        self.atom_number , self.atom_orb_number = States(self.projwfc_out_file)
        
        self.bandin_file_path = os.path.join(
            self.file_dir, self.bandin_file_name)
        self.bandin_file = open_file(self.bandin_file_path)
        self.nspin = n_spin(self.bandin_file)
        self.k_points_letter, self.k_points_bands, self.nk = k_point(self.bandin_file)
        
        

        self.verbosity = verbosity(self.bandin_file)
        if self.verbosity == 'high':
            self.bandout_file_path = os.path.join(self.file_dir, self.bandout_file_name)
            self.bandout_file = open_file(self.bandout_file_path)
            self.band_df= band_structure_from_bandsout(self.bandout_file, fermi_energy_corr, self.fermi_energy )

            points_break = [0]
            ct = 0
            for n in self.nk:
                ct +=n
                points_break.append(ct)
            k = []
            k_axis = self.band_df[self.band_df.columns[0]]
            self.k_axis = k_axis
            for i in range(len(self.band_df[self.band_df.columns[0]])):
                if (i in points_break) == True:
                    k.append(k_axis[i])
            self.k_points_path = array(k)
        else:
            print('verbosity must be high for projected bands.')
            
        
        pdos_files = []
        prefix = self.pdos_prefix
        
        for dir,subdir,files in os.walk(self.file_dir):
            for file in files:
                if ('pdos' in file) == True:
                    if file[:len(prefix)] == prefix and file != prefix+'.pdos.dat.pdos_tot' and file != prefix + '.pdos-proj.dat.projwfc_up' and file != prefix+'.pdos_tot':
                        pdos_files.append(os.path.join(self.file_dir,file))
                    elif file == prefix+'.pdos.dat.pdos_tot' or file == prefix+'.pdos_tot':
                        pdos_tot_file = os.path.join(self.file_dir,file)
                        
        def atm_wfc_number(file):
            for i in range(len(file)):
                if file[i] == 'm' and file[i+1] == '#':
                    n = i+3
                    atm = int(file[i+2])
                    break
            for i in range(n,len(file)):
                if file[i] == 'c' and file[i+1] == '#':
                    wfc = int(file[i+2])
                    break
            return [atm,wfc]
                    
        
        change_coluns_pdos(pdos_tot_file)
        
        pdos_tot = read_csv(pdos_tot_file, delim_whitespace= True,dtype=float)
        if fermi_energy_corr == True:
            pdos_tot[pdos_tot.columns[1]]-= self.fermi_energy
        pdos_tot.drop(pdos_tot.columns[3], axis='columns', inplace=True)

        
        self.pdos_kresolved = {0:pdos_tot}
        for file in pdos_files:
            change_coluns_pdos(file)
            atm,wfc = atm_wfc_number(file)
            if (atm in self.atom_orb_number.keys()) == True:
                if wfc == 1:
                    i = 3
                    for item in self.atom_orb_number[atm]['s']:
                        df = read_csv(file,delim_whitespace= True,dtype=float)
                        col = df.columns
                        if fermi_energy_corr == True:
                            df[col[1]]-=self.fermi_energy
                        for j in range(3,len(col)):
                            if j != i:
                                df.drop(col[j],axis='columns', inplace=True)
                        k = int(item)
                        df.drop(col[2],axis='columns', inplace=True)
                        self.pdos_kresolved[k] = df
                        i+=1
        
                elif wfc == 2:
                    i = 3
                    for item in self.atom_orb_number[atm]['p']:
                        df = read_csv(file,delim_whitespace= True,dtype=float)
                        col = df.columns
                        if fermi_energy_corr == True:
                            df[col[1]]-=self.fermi_energy
                        for j in range(3,len(col)):
                            if j != i:
                                df.drop(col[j],axis='columns', inplace=True)
                        k = int(item)
                        df.drop(col[2],axis='columns', inplace=True)
                        self.pdos_kresolved[k] = df
                        i+=1
        
                elif wfc == 3:
                    i = 3
                    for item in self.atom_orb_number[atm]['d']:
                        df = read_csv(file,delim_whitespace= True,dtype=float)
                        col = df.columns
                        if fermi_energy_corr == True:
                            df[col[1]]-=self.fermi_energy
                        for j in range(3,len(col)):
                            if j != i:
                                df.drop(col[j],axis='columns', inplace=True)
                        k = int(item)
                        df.drop(col[2],axis='columns', inplace=True)
                        self.pdos_kresolved[k] = df
                        i+=1
        
                elif wfc == 4:
                    i = 3
                    for item in self.atom_orb_number[atm]['f']:
                        df = read_csv(file,delim_whitespace= True,dtype=float)
                        col = df.columns
                        if fermi_energy_corr == True:
                            df[col[1]]-=self.fermi_energy
                        for j in range(3,len(col)):
                            if j != i:
                                df.drop(col[j],axis='columns', inplace=True)
                        k = int(item)
                        df.drop(col[2],axis='columns', inplace=True)
                        self.pdos_kresolved[k] = df
                        i+=1
        self.pdos_kresolved_plot = [self.pdos_kresolved,self.k_axis,self.k_points_path,self.k_points_letter]

    def orbitals_numbers(self):
        print(self.atom_number)
        print(self.atom_orb_number)

    def Charge2D(self):
        self.ppin_file_path = os.path.join(self.file_dir, self.ppin_file_name)
        self.ppin_file = open_file(self.ppin_file_path)
        self.parameters_pp = get_vectors(self.ppin_file)

        self.charge2D_file_path = os.path.join(self.file_dir, self.charge2D_name)
        
        def heatmat2D(df,nx,ny):
            X=[]
            Y=[]
            Z=[]
            for j in range(nx):
                teste_x=[]
                teste_y=[]
                teste_z=[]
                for i in range(ny):
                    teste_x.append(df.loc[i+j*nx][0])
                    teste_y.append(df.loc[i+j*nx][1])
                    teste_z.append(df.loc[i+j*nx][2])
                X.append(teste_x)
                Y.append(teste_y)
                Z.append(teste_z)
            return [X,Y,Z]
        
        df = heatmat2D(read_csv(self.charge2D_file_path, names= ['x','y','CD'],delim_whitespace=True),self.parameters_pp[3],self.parameters_pp[4])
        
        self.charge2D_data = [df,self.alat]
        


class Plot:
    def bands(self, data, k_points_letter, k_path, subplotsize=(10, 8), subplot=False, ax=None, vline=True, vline_linewidth=1, E_min=-5, E_max=5, dE=1, font_axis=20, font_label=20, ylabel='E -EF (eV)', xlabel='', title='Band Structure',
              linewidth=3, legend=False, loc="upper right", color='black', color_vline='gray', vline_linestyle='dashed', label_band='eletronic strucuture', fermi_level_line=False, fermi_level=0, fermi_level_color='red', 
              fermi_level_linewidth=1, fermi_level_linestyle='dashed', occupied_and_unoccupied_states_dot = False, occupied_and_unoccupied_states_points = ((0,0),(0,0)), scatter_type = 'o',occupied_and_unoccupied_states_color_dot = 'blue',
              occupied_and_unoccupied_states_line=False, occupied_and_unoccupied_states=(-1, 1), occupied_and_unoccupied_states_color=('blue', 'blue'), occupied_and_unoccupied_states_linewidth=1,occupied_and_unoccupied_states_size_dot=80
              , dpi = 150) :
        # flag subplot
        if subplot == False:
            del (ax)
            fig, ax = plt.subplots(figsize=subplotsize, dpi = dpi)
        # -----------------------------------------------------------------------------------------------------
        # plot data
        columns = data.columns
        k = data[columns[0]]
        for i in range(1, len(columns)):
            energy = data[columns[i]]
            if legend == True:
                if i == 1:
                    ax.plot(k, energy, color=color,
                            linewidth=linewidth, label=label_band)
                else:
                    ax.plot(k, energy, color=color, linewidth=linewidth)
            else:
                ax.plot(k, energy, color=color, linewidth=linewidth)
        # -----------------------------------------------------------------------------------------------------
        # flag vline
        if vline == True:
            for i in range(1, len(k_path)-1):
                ax.vlines(k_path[i], E_min, E_max, linestyles=vline_linestyle,
                          colors=color_vline, linewidth=vline_linewidth)
        # -----------------------------------------------------------------------------------------------------
        # flag fermy_level_line
        if fermi_level_line == True:
            ax.hlines(fermi_level, k_path[0], k_path[len(
                k_path)-1], linestyles=fermi_level_linestyle, color=fermi_level_color, linewidth=fermi_level_linewidth)
        # -----------------------------------------------------------------------------------------------------
        # flag occupied_and_unoccupied_states_line
        if occupied_and_unoccupied_states_line == True and occupied_and_unoccupied_states_dot == False:
            for i in range(len(occupied_and_unoccupied_states)):
                item = occupied_and_unoccupied_states[i]
                ax.hlines(item, k_path[0], k_path[len(k_path)-1], linestyles=fermi_level_linestyle,
                          linewidth=occupied_and_unoccupied_states_linewidth, color=occupied_and_unoccupied_states_color[i])
        elif occupied_and_unoccupied_states_line == False and occupied_and_unoccupied_states_dot == True:
            x = []; y = []
            for item in occupied_and_unoccupied_states_points:
                x.append(item[0]); y.append(item[1])
            ax.scatter(x,y, marker = scatter_type,s= occupied_and_unoccupied_states_size_dot, color = occupied_and_unoccupied_states_color_dot)
        # -----------------------------------------------------------------------------------------------------

        # flag legend
        if legend == True:
            ax.legend(loc=loc, fontsize=font_label)
        # -----------------------------------------------------------------------------------------------------
        # globar parameters
        ax.set_xlim(k_path[0], k_path[len(k_path)-1])
        ax.set_xticks(k_path, k_points_letter, fontsize=font_axis)
        ax.set_xlabel(xlabel)

        ax.set_ylim(E_min, E_max)
        ax.set_yticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)
        ax.set_ylabel(ylabel)
        ax.xaxis.label.set_size(font_label)
        ax.yaxis.label.set_size(font_label)

        ax.set_title(title, fontsize=font_label)

        return ax

    def dos(self, data, subplotsize=(10, 8), subplot=False, ax=None, color='black', linewidth=3, label='dos', legend = False, loc = 'upper right', font_label = 20, fill = False, fill_color ='blue', alpha = 0.5,
            E_min = -5, E_max = 5 ,dE = 1, dos_min = 0, dos_max = 10, delta_dos = 2, font_axis = 20, xlabel = 'E - E_F (eV)', ylabel = 'DoS (states/eV)', title = 'Density of States',
            fermi_level_line = False, fermi_level= 0, fermi_level_color='red', fermi_level_linewidth=1, fermi_level_linestyle='dashed'
            ):
        # flag subplot
        if subplot == False:
            del (ax)
            fig, ax = plt.subplots(figsize=subplotsize)
        # plot data
        energy = data[data.columns[0]]
        dos = data[data.columns[1]]
        ax.plot(energy, dos, color = color, linewidth = linewidth, label = label)

        #flag fill_between
        if fill == True:
            ax.fill_between(energy,dos, color = fill_color   ,alpha = alpha)

        # flag fermy_level_line
        if fermi_level_line == True:
            ax.vlines(fermi_level, dos_min, dos_max , linestyles=fermi_level_linestyle, color=fermi_level_color, linewidth=fermi_level_linewidth)

        # flag legend
        if legend == True:
            ax.legend(loc=loc, fontsize=font_label)

        # globar parameters
        ax.set_xlim(E_min, E_max)
        ax.set_xticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)
        ax.set_xlabel(xlabel)

        ax.set_ylim(dos_min, dos_max)
        ax.set_yticks(arange(dos_min, dos_max+delta_dos, delta_dos),
                      arange(dos_min, dos_max+delta_dos, delta_dos), fontsize=font_axis)
        ax.set_ylabel(ylabel)
        ax.xaxis.label.set_size(font_label)
        ax.yaxis.label.set_size(font_label)

        ax.set_title(title, fontsize=font_label)
        return ax

    def pdos_atoms(self, dict, subplotsize=(10, 8), subplot=False, ax=None, all_atom = True, atom = ['A','B'], linestyle = 'solid',
     linewidth=3, legend = True, loc = 'upper right', font_label = 20, fill = False, fill_color ='blue', alpha = 0.5,
            E_min = -5, E_max = 5 ,dE = 1, dos_min = 0, dos_max = 10, delta_dos = 1, font_axis = 20, xlabel = 'E - E_F (eV)', ylabel = 'DoS (states/eV)', title = 'Density of States per Atom',
            fermi_level_line = False, fermi_level= 0, fermi_level_color='red', fermi_level_linewidth=1, fermi_level_linestyle='dashed'):
        #flag subplot
        if subplot == False:
            del (ax)
            fig, ax = plt.subplots(figsize=subplotsize)
    
        list_data = []
        for key in dict.keys():
            ### flag all_atom
            if all_atom == True:
                data = dict[key]
                x = data[data.columns[0]]
                y = data[data.columns[1]]
                ax.plot(x, y, linewidth = linewidth, label = key , linestyle = linestyle)
                list_data.append(y)
            elif all_atom == False:
                if (key in atom) == True:
                    data = dict[key] 
                    x = data[data.columns[0]]
                    y = data[data.columns[1]]
                    list_data.append(y)
                    ax.plot(x, y, linewidth = linewidth, label = key ,linestyle = linestyle)
            else:
                print('Wrong value for all_atom parameter, only accepet True or False. Default = True.')

        
        #flag fill_between
        #if fill == True:
        #    ax.fill_between(energy,dos, color = fill_color   ,alpha = alpha)

        # flag fermy_level_line
        if fermi_level_line == True:
            ax.vlines(fermi_level, dos_min, dos_max , linestyles=fermi_level_linestyle, color=fermi_level_color, linewidth=fermi_level_linewidth)

        # flag legend
        if legend == True:
            ax.legend(loc=loc, fontsize=font_label)

        # globar parameters
        ax.set_xlim(E_min, E_max)
        ax.set_xticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)
        ax.set_xlabel(xlabel)

        ax.set_ylim(dos_min, dos_max)
        ax.set_yticks(arange(dos_min, dos_max+delta_dos, delta_dos),
                      arange(dos_min, dos_max+delta_dos, delta_dos), fontsize=font_axis)
        ax.set_ylabel(ylabel)
        ax.xaxis.label.set_size(font_label)
        ax.yaxis.label.set_size(font_label)

        ax.set_title(title, fontsize=font_label)
        return ax

    def pdos_orb(self, dict, subplotsize=(10, 8), subplot=False, ax=None, all_orb = True, orb = {'A':'all','B':'all'},  linestyle = 'solid',
     linewidth=3, legend = True, loc = 'upper right', font_label = 20, fill = False, fill_color ='blue', alpha = 0.5,
            E_min = -5, E_max = 5 ,dE = 1, dos_min = 0, dos_max = 10, delta_dos = 1, font_axis = 20, xlabel = 'E - E_F (eV)', ylabel = 'DoS (states/eV)', title = 'Density of States per Atom',
            fermi_level_line = False, fermi_level= 0, fermi_level_color='red', fermi_level_linewidth=1, fermi_level_linestyle='dashed'):
        if subplot == False:
            del (ax)
            fig, ax = plt.subplots(figsize=subplotsize)
        for key1 in dict.keys():
            for key2 in dict[key1].keys():
                if all_orb == True:
                    data = dict[key1][key2]
                    x = data[data.columns[0]]
                    y = data[data.columns[1]]
                    ax.plot(x, y, label=key1+' ' +key2 , linewidth = linewidth, linestyle = linestyle)
                elif all_orb == False:
                    if (key1 in orb.keys()) == True:
                        if orb[key1] == 'all':
                            data = dict[key1][key2]
                            x = data[data.columns[0]]
                            y = data[data.columns[1]]
                            ax.plot(x, y, label=key1+' ' +key2,linewidth = linewidth,linestyle = linestyle)
                        elif (key2 in orb[key1]) == True:
                            data = dict[key1][key2]
                            x = data[data.columns[0]]
                            y = data[data.columns[1]]
                            ax.plot(x, y, label=key1+ ' ' + key2,linewidth = linewidth,linestyle = linestyle)                            

                else:
                    print('Wrong value for all_orb parameter, only accepet True or False. Default = True.')

                # flag fermy_level_line
        if fermi_level_line == True:
            ax.vlines(fermi_level, dos_min, dos_max , linestyles=fermi_level_linestyle, color=fermi_level_color, linewidth=fermi_level_linewidth)

        # flag legend
        if legend == True:
            ax.legend(loc=loc, fontsize=font_label)

        # globar parameters
        ax.set_xlim(E_min, E_max)
        ax.set_xticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)
        ax.set_xlabel(xlabel)

        ax.set_ylim(dos_min, dos_max)
        ax.set_yticks(arange(dos_min, dos_max+delta_dos, delta_dos),
                      arange(dos_min, dos_max+delta_dos, delta_dos), fontsize=font_axis)
        ax.set_ylabel(ylabel)
        ax.xaxis.label.set_size(font_label)
        ax.yaxis.label.set_size(font_label)

        ax.set_title(title, fontsize=font_label)
        return ax

    def bands_dos(self, data_band, k_points_letter, k_path, data_dos, subplotsize=[10, 8], vline_bands=True, vline_linewidth_bands=1, E_min=-5, E_max=5, dE=1, font_axis=20, font_label=20, ylabel_bands='E -EF (eV)', xlabel_bands='', title_bands='Band Structure',
              linewidth_bands=3, legend_bands=False, loc_bands="upper right", color_bands='black', color_vline_bands='gray', vline_linestyle_bands='dashed', label_band_bands='eletronic strucuture', fermi_level_line_bands=False, fermi_level=0, fermi_level_color_bands='red',
              fermi_level_linewidth_bands=1, fermi_level_linestyle_bands='dashed',
              occupied_and_unoccupied_states_line=False, occupied_and_unoccupied_states=[-1, 1], occupied_and_unoccupied_states_color=('blue', 'blue'), occupied_and_unoccupied_states_linewidth=1,
              occupied_and_unoccupied_states_dot = False, occupied_and_unoccupied_states_points = ((0,0),(0,0)), scatter_type = 'o',occupied_and_unoccupied_states_color_dot = 'blue', occupied_and_unoccupied_states_size_dot=80 ,        
              
            linewidth_dos=3, label_dos='dos', legend_dos = True, loc_dos = 'upper right', font_label_dos = 20, fill = True, fill_color ='blue', alpha = 0.5, color_dos = 'black',
            dos_max = 10, delta_dos = 2, font_axis_dos = 20, xlabel_dos = 'DoS (states/eV)', title_dos = 'Density of States',
            fermi_level_line_dos = True, fermi_level_dos= 0, fermi_level_color_dos='red', fermi_level_linewidth_dos=1, fermi_level_linestyle_dos='dashed',

              left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0, hspace=0.1
              ):
                
        fig, ax = plt.subplots(1,2,figsize=subplotsize)

        # -----------------------------------------------------------------------------------------------------
        # plot data
        columns = data_band.columns
        k = data_band[columns[0]]
        for i in range(1, len(columns)):
            energy = data_band[columns[i]]
            if legend_bands == True:
                if i == 1:
                    ax[0].plot(k, energy, color=color_bands,
                            linewidth=linewidth_bands, label=label_band_bands)
                else:
                    ax[0].plot(k, energy, color=color_bands, linewidth=linewidth_bands)
            else:
                ax[0].plot(k, energy, color=color_bands, linewidth=linewidth_bands)
        # -----------------------------------------------------------------------------------------------------
        # flag vline
        if vline_bands == True:
            for i in range(1, len(k_path)-1):
                ax[0].vlines(k_path[i], E_min, E_max, linestyles=vline_linestyle_bands,
                          colors=color_vline_bands, linewidth=vline_linewidth_bands)
        # -----------------------------------------------------------------------------------------------------
        # flag fermy_level_line
        if fermi_level_line_bands == True:
            ax[0].hlines(fermi_level, k_path[0], k_path[len(
                k_path)-1], linestyles=fermi_level_linestyle_bands, color=fermi_level_color_bands, linewidth=fermi_level_linewidth_bands)
        # -----------------------------------------------------------------------------------------------------
        if occupied_and_unoccupied_states_line == True and occupied_and_unoccupied_states_dot == False:
            for i in range(len(occupied_and_unoccupied_states)):
                item = occupied_and_unoccupied_states[i]
                ax[0].hlines(item, k_path[0], k_path[len(k_path)-1], linestyles=fermi_level_linestyle_bands,
                          linewidth=occupied_and_unoccupied_states_linewidth, color=occupied_and_unoccupied_states_color[i])
        elif occupied_and_unoccupied_states_line == False and occupied_and_unoccupied_states_dot == True:
            x = []; y = []
            for item in occupied_and_unoccupied_states_points:
                x.append(item[0]); y.append(item[1])
            ax[0].scatter(x,y, marker = scatter_type,s = occupied_and_unoccupied_states_size_dot, color = occupied_and_unoccupied_states_color_dot)
        # -----------------------------------------------------------------------------------------------------

        # flag legend
        if legend_bands == True:
            ax[0].legend(loc=loc_bands, fontsize=font_label)
        # -----------------------------------------------------------------------------------------------------
        # globar parameters
        ax[0].set_xlim(k_path[0], k_path[len(k_path)-1])
        ax[0].set_xticks(k_path, k_points_letter, fontsize=font_axis)
        ax[0].set_xlabel(xlabel_bands)

        ax[0].set_ylim(E_min, E_max)
        ax[0].set_yticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)
        ax[0].set_ylabel(ylabel_bands)
        ax[0].xaxis.label.set_size(font_label)
        ax[0].yaxis.label.set_size(font_label)

        ax[0].set_title(title_bands, fontsize=font_label)

        #########################################################################################################
        #### plot dos

        # plot data
        energy = data_dos[data_dos.columns[0]]
        dos = data_dos[data_dos.columns[1]]
        ax[1].plot(dos, energy, color = color_dos, linewidth = linewidth_dos, label = label_dos)

        #flag fill_between
        if fill == True:
            ax[1].fill_betweenx(energy,dos ,color = fill_color   ,alpha = alpha)

        # flag fermy_level_line
        if fermi_level_line_dos == True:
            ax[1].hlines(fermi_level, 0, dos_max , linestyles=fermi_level_linestyle_dos, color=fermi_level_color_dos, linewidth=fermi_level_linewidth_dos)

        # flag legend
        if legend_dos == True:
            ax[1].legend(loc=loc_dos, fontsize=font_label)

        # globar parameters
        ax[1].set_ylim(E_min, E_max)
        ax[1].set_yticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)

        ax[1].set_xlim(0, dos_max)
        ax[1].set_xticks(arange(delta_dos, dos_max+delta_dos, delta_dos),
                      arange(delta_dos, dos_max+delta_dos, delta_dos), fontsize=font_axis)
        ax[1].set_xlabel(xlabel_dos)
        ax[1].xaxis.label.set_size(font_label)
        ax[1].yaxis.label.set_size(font_label)

        ax[1].set_title(title_dos, fontsize=font_label)

        ax[1].yaxis.tick_right()
    

        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        return ax

    def bands_pdos_atom(self, data_band, k_points_letter, k_path, dict, subplotsize=[10, 8], vline_bands=True, vline_linewidth_bands=1, E_min=-5, E_max=5, dE=1, font_axis=20, font_label=20, ylabel_bands='E -EF (eV)', xlabel_bands='', title_bands='Band Structure',
              linewidth_bands=3, legend_bands=False, loc_bands="upper right", color_bands='black', color_vline_bands='gray', vline_linestyle_bands='dashed', label_band_bands='eletronic strucuture', fermi_level_line_bands=False, fermi_level=0, fermi_level_color_bands='red',
              fermi_level_linewidth_bands=1, fermi_level_linestyle_bands='dashed',
              occupied_and_unoccupied_states_line=False, occupied_and_unoccupied_states=[-1, 1], occupied_and_unoccupied_states_color=('blue', 'blue'), occupied_and_unoccupied_states_size_dot=80,occupied_and_unoccupied_states_linewidth =1,
              occupied_and_unoccupied_states_dot = False, occupied_and_unoccupied_states_points = ((0,0),(0,0)), scatter_type = 'o',occupied_and_unoccupied_states_color_dot = 'blue',          
            all_atom = True, atom = ['A','B'], linestyle_pdos = 'solid',
        linewidth_pdos=3, legend_pdos = True, loc_pdos = 'upper right',
            dos_max = 10, delta_dos = 2, xlabel_pdos = 'DoS (states/eV)',title_pdos = 'Density of States per Atom',
            fermi_level_line_pdos = False, fermi_level_color_pdos='red', fermi_level_linewidth_pdos=1, fermi_level_linestyle_pdos='dashed',

              left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0, hspace=0.1
              ):
                
        fig, ax = plt.subplots(1,2,figsize=subplotsize)

        # -----------------------------------------------------------------------------------------------------
        # plot data
        columns = data_band.columns
        k = data_band[columns[0]]
        for i in range(1, len(columns)):
            energy = data_band[columns[i]]
            if legend_bands == True:
                if i == 1:
                    ax[0].plot(k, energy, color=color_bands,
                            linewidth=linewidth_bands, label=label_band_bands)
                else:
                    ax[0].plot(k, energy, color=color_bands, linewidth=linewidth_bands)
            else:
                ax[0].plot(k, energy, color=color_bands, linewidth=linewidth_bands)
        # -----------------------------------------------------------------------------------------------------
        # flag vline
        if vline_bands == True:
            for i in range(1, len(k_path)-1):
                ax[0].vlines(k_path[i], E_min, E_max, linestyles=vline_linestyle_bands,
                          colors=color_vline_bands, linewidth=vline_linewidth_bands)
        # -----------------------------------------------------------------------------------------------------
        # flag fermy_level_line
        if fermi_level_line_bands == True:
            ax[0].hlines(fermi_level, k_path[0], k_path[len(
                k_path)-1], linestyles=fermi_level_linestyle_bands, color=fermi_level_color_bands, linewidth=fermi_level_linewidth_bands)
        # -----------------------------------------------------------------------------------------------------
        if occupied_and_unoccupied_states_line == True and occupied_and_unoccupied_states_dot == False:
            for i in range(len(occupied_and_unoccupied_states)):
                item = occupied_and_unoccupied_states[i]
                ax[0].hlines(item, k_path[0], k_path[len(k_path)-1], linestyles=fermi_level_linestyle_bands,
                          linewidth=occupied_and_unoccupied_states_linewidth, color=occupied_and_unoccupied_states_color[i])
        elif occupied_and_unoccupied_states_line == False and occupied_and_unoccupied_states_dot == True:
            x = []; y = []
            for item in occupied_and_unoccupied_states_points:
                x.append(item[0]); y.append(item[1])
            ax[0].scatter(x,y, marker = scatter_type,s =  occupied_and_unoccupied_states_size_dot, color = occupied_and_unoccupied_states_color_dot)
        # -----------------------------------------------------------------------------------------------------

        # flag legend
        if legend_bands == True:
            ax[0].legend(loc=loc_bands, fontsize=font_label)
        # -----------------------------------------------------------------------------------------------------
        # globar parameters
        ax[0].set_xlim(k_path[0], k_path[len(k_path)-1])
        ax[0].set_xticks(k_path, k_points_letter, fontsize=font_axis)
        ax[0].set_xlabel(xlabel_bands)

        ax[0].set_ylim(E_min, E_max)
        ax[0].set_yticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)
        ax[0].set_ylabel(ylabel_bands)
        ax[0].xaxis.label.set_size(font_label)
        ax[0].yaxis.label.set_size(font_label)

        ax[0].set_title(title_bands, fontsize=font_label)

        #########################################################################################################
        #### plot pdos

        # plot data
        list_data = []
        for key in dict.keys():
            ### flag all_atom
            if all_atom == True:
                data = dict[key]
                x = data[data.columns[0]]
                y = data[data.columns[1]]
                ax[1].plot(y, x, linewidth = linewidth_pdos, label = key , linestyle = linestyle_pdos)
                list_data.append(y)
            elif all_atom == False:
                if  (key in atom) == True:
                    data = dict[key] 
                    x = data[data.columns[0]]
                    y = data[data.columns[1]]
                    list_data.append(y)
                    ax[1].plot(y, x, linewidth = linewidth_pdos, label = key ,linestyle = linestyle_pdos)
            else:
                print('Wrong value for all_atom parameter, only accepet True or False. Default = True.')

        
        #flag fill_between
        #if fill == True:
        #    ax.fill_between(energy,dos, color = fill_color   ,alpha = alpha)

        # flag fermy_level_line
        if fermi_level_line_pdos == True:
            ax[1].vlines(fermi_level,0, dos_max , linestyles=fermi_level_linestyle_pdos, color=fermi_level_color_pdos, linewidth=fermi_level_linewidth_pdos)

        # flag legend
        if legend_pdos == True:
            ax[1].legend(loc=loc_pdos, fontsize=font_label)

        # globar parameters
        ax[1].set_ylim(E_min, E_max)
        ax[1].set_yticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)


        ax[1].set_xlim(0, dos_max)
        ax[1].set_xticks(arange(delta_dos, dos_max+delta_dos, delta_dos),
                      arange(delta_dos, dos_max+delta_dos, delta_dos), fontsize=font_axis)

        ax[1].set_xlabel(xlabel_pdos)
        ax[1].xaxis.label.set_size(font_label)
        ax[1].yaxis.label.set_size(font_label)

        ax[1].set_title(title_pdos, fontsize=font_label)

        ax[1].yaxis.tick_right()
    

        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        return ax

    def bands_pdos_orb(self, data_band, k_points_letter, k_path, dict, subplotsize=[10, 8],dpi = 150, vline_bands=True, vline_linewidth_bands=1, E_min=-5, E_max=5, dE=1, font_axis=20, font_label=20, ylabel_bands='E -EF (eV)', xlabel_bands='', title_bands='Band Structure',
              linewidth_bands=3, legend_bands=False, loc_bands="upper right", color_bands='black', color_vline_bands='gray', vline_linestyle_bands='dashed', label_band_bands='eletronic strucuture', fermi_level_line_bands=False, fermi_level=0, fermi_level_color_bands='red',
              fermi_level_linewidth_bands=1, fermi_level_linestyle_bands='dashed',
              occupied_and_unoccupied_states_line=False, occupied_and_unoccupied_states=[-1, 1], occupied_and_unoccupied_states_color=('blue', 'blue'), occupied_and_unoccupied_states_linewidth=1,
              occupied_and_unoccupied_states_dot = False, occupied_and_unoccupied_states_points = ((0,0),(0,0)), scatter_type = 'o',occupied_and_unoccupied_states_color_dot = 'blue',          
            all_orb = True, orb = ['A','B'], linestyle_pdos = 'solid',
        linewidth_pdos=3, legend_pdos = True, loc_pdos = 'upper right',
            dos_max = 10, delta_dos = 2, xlabel_pdos = 'DoS (states/eV)',title_pdos = 'Density of States per Atom',
            fermi_level_line_pdos = False, fermi_level_color_pdos='red', fermi_level_linewidth_pdos=1, fermi_level_linestyle_pdos='dashed',

              left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0, hspace=0.1
              ):
                
        fig, ax = plt.subplots(1,2,figsize=subplotsize, dpi = dpi)

        # -----------------------------------------------------------------------------------------------------
        # plot data
        columns = data_band.columns
        k = data_band[columns[0]]
        for i in range(1, len(columns)):
            energy = data_band[columns[i]]
            if legend_bands == True:
                if i == 1:
                    ax[0].plot(k, energy, color=color_bands,
                            linewidth=linewidth_bands, label=label_band_bands)
                else:
                    ax[0].plot(k, energy, color=color_bands, linewidth=linewidth_bands)
            else:
                ax[0].plot(k, energy, color=color_bands, linewidth=linewidth_bands)
        # -----------------------------------------------------------------------------------------------------
        # flag vline
        if vline_bands == True:
            for i in range(1, len(k_path)-1):
                ax[0].vlines(k_path[i], E_min, E_max, linestyles=vline_linestyle_bands,
                          colors=color_vline_bands, linewidth=vline_linewidth_bands)
        # -----------------------------------------------------------------------------------------------------
        # flag fermy_level_line
        if fermi_level_line_bands == True:
            ax[0].hlines(fermi_level, k_path[0], k_path[len(
                k_path)-1], linestyles=fermi_level_linestyle_bands, color=fermi_level_color_bands, linewidth=fermi_level_linewidth_bands)
        # -----------------------------------------------------------------------------------------------------
        if occupied_and_unoccupied_states_line == True and occupied_and_unoccupied_states_dot == False:
            for i in range(len(occupied_and_unoccupied_states)):
                item = occupied_and_unoccupied_states[i]
                ax[0].hlines(item, k_path[0], k_path[len(k_path)-1], linestyles=fermi_level_linestyle_bands,
                          linewidth=occupied_and_unoccupied_states_linewidth, color=occupied_and_unoccupied_states_color[i])
        elif occupied_and_unoccupied_states_line == False and occupied_and_unoccupied_states_dot == True:
            x = []; y = []
            for item in occupied_and_unoccupied_states_points:
                x.append(item[0]); y.append(item[1])
            ax[0].plot(x,y, scatter_type,linewidth=occupied_and_unoccupied_states_linewidth, color = occupied_and_unoccupied_states_color_dot)
        # -----------------------------------------------------------------------------------------------------

        # flag legend
        if legend_bands == True:
            ax[0].legend(loc=loc_bands, fontsize=font_label)
        # -----------------------------------------------------------------------------------------------------
        # globar parameters
        ax[0].set_xlim(k_path[0], k_path[len(k_path)-1])
        ax[0].set_xticks(k_path, k_points_letter, fontsize=font_axis)
        ax[0].set_xlabel(xlabel_bands)

        ax[0].set_ylim(E_min, E_max)
        ax[0].set_yticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)
        ax[0].set_ylabel(ylabel_bands)
        ax[0].xaxis.label.set_size(font_label)
        ax[0].yaxis.label.set_size(font_label)

        ax[0].set_title(title_bands, fontsize=font_label)

        #########################################################################################################
        #### plot pdos

        # plot data
        for key1 in dict.keys():
            for key2 in dict[key1].keys():
                if all_orb == True:
                    data = dict[key1][key2]
                    x = data[data.columns[0]]
                    y = data[data.columns[1]]
                    ax[1].plot(y, x, label=key1+' ' +key2)
                elif all_orb == False:
                    if (key1 in orb.keys()) == True:
                        if orb[key1] == 'all':
                            data = dict[key1][key2]
                            x = data[data.columns[0]]
                            y = data[data.columns[1]]
                            ax[1].plot(y, x, label=key1+' ' +key2)
                        elif (key2 in orb[key1]) == True:
                            data = dict[key1][key2]
                            x = data[data.columns[0]]
                            y = data[data.columns[1]]
                            ax[1].plot(y, x, label=key1+ ' ' + key2)                            

                else:
                    print('Wrong value for all_orb parameter, only accepet True or False. Default = True.')

        
        #flag fill_between
        #if fill == True:
        #    ax.fill_between(energy,dos, color = fill_color   ,alpha = alpha)

        # flag fermy_level_line
        if fermi_level_line_pdos == True:
            ax[1].vlines(fermi_level,0, dos_max , linestyles=fermi_level_linestyle_pdos, color=fermi_level_color_pdos, linewidth=fermi_level_linewidth_pdos)

        # flag legend
        if legend_pdos == True:
            ax[1].legend(loc=loc_pdos, fontsize=font_label)

        # globar parameters
        ax[1].set_ylim(E_min, E_max)
        ax[1].set_yticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)


        ax[1].set_xlim(0, dos_max)
        ax[1].set_xticks(arange(delta_dos, dos_max+delta_dos, delta_dos),
                      arange(delta_dos, dos_max+delta_dos, delta_dos), fontsize=font_axis)

        ax[1].set_xlabel(xlabel_pdos)
        ax[1].xaxis.label.set_size(font_label)
        ax[1].yaxis.label.set_size(font_label)

        ax[1].set_title(title_pdos, fontsize=font_label)

        ax[1].yaxis.tick_right()
    

        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        return ax

    def bands_project(self,bands_proj, subplotsize=(10, 8), subplot = False, ax = None,dpi = 150,
    legend = False, font_label = 20, font_axis = 20, xlabel = '', ylabel = 'E - E_f (eV)', E_min = -5, E_max = 5, dE = 1,loc = 'upper right',title = 'Band projected',
    vline = True, color_vline='gray', vline_linestyle='dashed', vline_linewidth = 1.5, backgorund = 'white', cmap = 'viridis', norm = (0,1)
    , index = False, states = [0,1,2], corretalion = False, states1 = [1,2,3], states2 = [4,5,6],norm_corr = (-1,1), cbar_label = ['0','1'], cbar_tick = [0,1]
    ,linewidth = 3):


        if subplot == False:
            del (ax)
            fig, ax = plt.subplots(figsize=subplotsize, dpi = dpi)
        bands = bands_proj[0]; k_letter = bands_proj[2]; k_path = bands_proj[3]
        columns = bands.columns
        x_interp = bands[columns[0]]
        x = arange(0,x_interp.max(),0.001)
        for i in range(1,len(columns)):
            y_int = bands[columns[i]]
            if corretalion == False:
                proj = bands_proj[1]
                marker = False
                for state in states:
                    if marker == False:
                        z_int = array(proj[state][i])
                        marker = True
                    else:
                        z_int= z_int  + array(proj[state][i])
                        
            elif corretalion == True:
                proj = bands_proj[1]
                marker = False
                for state1 in states1:
                    if marker == False:
                        z1 = array(proj[state1][i])
                        marker = True
                    else:
                        z1= z1  + array(proj[state1][i])
                marker = False
                for state2 in states2:
                    if marker == False:
                        z2 = array(proj[state1][i])
                        marker = True
                    else:
                        z2= z2  + array(proj[state1][i])
                z_int = z2 - z1
            else:
                print('Wrong value for correlation, only True and False are accepted. Default is False')
            
                
            f = interpolate.interp1d(x_interp, y_int)
            g = interpolate.interp1d(x_interp, z_int)
            y = f(x)
            z = g(x)
            points = array([x, y]).T.reshape(-1, 1, 2)
            segments = concatenate([points[:-1], points[1:]], axis=1)
            if corretalion == False:
                norma = plt.Normalize(norm[0],norm[1])
            elif corretalion == True:
                norma = plt.Normalize(norm_corr[0],norm_corr[1])
            #norm = plt.Normalize(0, z.max())
            #norm = plt.Normalize(mini,maxi)
            lc = LineCollection(segments, cmap=cmap, norm=norma, linewidth = linewidth) #cmap='viridis'
            lc.set_array(z)
            lc.set_linewidth(2)
            line = ax.add_collection(lc)

        ax.set_facecolor(backgorund)
        cbar = fig.colorbar(line, ticks = cbar_tick)
        cbar.ax.set_yticklabels(cbar_label,fontsize = 20)


        # flag vline
        if vline == True:
            for i in range(1, len(k_path)-1):
                ax.vlines(k_path[i], E_min, E_max, linestyles=vline_linestyle,
                          colors=color_vline, linewidth=vline_linewidth)
        # -----------------------------------------------------------------------------------------------------

                # flag legend
        if legend == True:
            ax.legend(loc=loc, fontsize=font_label)
        # -----------------------------------------------------------------------------------------------------
        # globar parameters
        ax.set_xlim(k_path[0], k_path[len(k_path)-1])
        ax.set_xticks(k_path, k_letter, fontsize=font_axis)
        ax.set_xlabel(xlabel)

        ax.set_ylim(E_min, E_max)
        ax.set_yticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)
        ax.set_ylabel(ylabel)
        ax.xaxis.label.set_size(font_label)
        ax.yaxis.label.set_size(font_label)

        ax.set_title(title, fontsize=font_label)

        return ax
        
    def kresolved(self, data, states = [0],subplotsize=(10, 8), subplot=False, ax=None, ylabel = 'E - E_f (eV)', E_min = -5, E_max = 5, dE = 1,loc = 'upper right',title = 'Band projected'
                  , cmap = 'magma', font_axis = 20, font_label = 20,cbar_label = ['0','1'], cbar_tick = [0,1]):
        
        pdos_kresolved = data[0];k_axis = data[1]; k_path = data[2]; k_letter =data[3]
        if subplot == False:
            del (ax)
            fig, ax = plt.subplots(figsize=subplotsize)
            
        marker= False
        for state in states:
            if marker == False:
                df = pdos_kresolved[state]
            
                column = df.columns
                df.rename(columns={column[2]: 'ldos'})
                column = df.columns
                k_tot = df[column[0]]; e_tot = df[column[1]] 
                dos_tot = df[column[2]]
                marker = True
            else:
                df = pdos_kresolved[state]
                column = df.columns
                df.rename(columns={column[2]: 'ldos'})
                column = df.columns
                dos_tot = dos_tot + df[column[2]]
            
        n = 0
        e = []
        for i in range(len(k_tot)):
            if k_tot[i] == 1:
                e.append(e_tot[i])
                n+=1
                
        dos = []; K = []; E = []
        for i in range(len(k_axis)):
            col = []
            for j in range(n):
                col.append(k_axis[i])
            K.append(col)
            E.append(e_tot[i*n:(i+1)*n])
            dos.append(array(dos_tot[i*n:(i+1)*n])/dos_tot.max())

        im = ax.pcolormesh(K,E,dos, cmap=cmap, shading='auto')


        # globar parameters
        ax.set_xlim(k_path[0], k_path[len(k_path)-1])
        ax.set_xticks(k_path, k_letter, fontsize=font_axis)
        
        ax.set_ylim(E_min, E_max)
        ax.set_yticks(arange(E_min, E_max+dE, dE),
                      arange(E_min, E_max+dE, dE), fontsize=font_axis)
        ax.set_ylabel(ylabel)
        ax.xaxis.label.set_size(font_label)
        ax.yaxis.label.set_size(font_label)

        ax.set_title(title, fontsize=font_label)
        
        cbar = fig.colorbar(im, ticks = cbar_tick)
        cbar.ax.set_yticklabels(cbar_label,fontsize = font_axis)
        
    def charge_density2D(self,data,subplotsize=(10, 8), subplot=False, ax=None,xlim = (0,10,5),ylim = (0,10,5), font_axis = 20, cbar_label = ['0','+'], cbar_tick = [0,0.5]):
        
        chargedensity2D = data[0]; alat_bohr = data[1];alat_ang = data[1]/1.88973
        if subplot == False:
            del (ax)
            fig, ax = plt.subplots(figsize=subplotsize)
        
        
        im = ax.pcolormesh(chargedensity2D[0],chargedensity2D[1] ,chargedensity2D[2])
        
 
        ax.set_xticks( arange( xlim[0]*1.88973, xlim[1]*1.88973 +  xlim[2]*1.88973, xlim[2]*1.88973),  arange( xlim[0], xlim[1] +  xlim[2], xlim[2]),fontsize=font_axis )
    
        
        ax.set_xlim(xlim[0]*1.88973,xlim[1]*1.88973 )
        
        
        
        ax.set_yticks( arange( ylim[0]*1.88973, ylim[1]*1.88973 +  ylim[2]*1.88973, ylim[2]*1.88973),  arange( ylim[0], ylim[1] +  ylim[2], ylim[2]), fontsize=font_axis)
    
        
        ax.set_ylim(ylim[0]*1.88973,ylim[1] *1.88973,)
        
        cbar = fig.colorbar(im, ticks = cbar_tick)
        cbar.ax.set_yticklabels(cbar_label,fontsize = font_axis)