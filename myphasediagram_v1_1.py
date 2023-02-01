from sympy import Symbol, nsolve, solve
from scipy.optimize import fsolve
import numpy as np
import math
import itertools
import pymatgen.core as mg
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.ext.matproj import MPRester
from sympy import Matrix
from fractions import Fraction, gcd
from functools import reduce
import plotly.graph_objs as go
import ternary
from matminer.data_retrieval.retrieve_MP import MPDataRetrieval
import pandas as pd

def get_phasediagram (key, elements):
    mpr = MPRester(key)
    entries = mpr.get_entries_in_chemsys(elements)
    #print(entries)
    pd = PhaseDiagram(entries)
    return pd



def get_specific_compounds (pd, e_above_hull):

    list_of_all = list()
    for ent in pd.all_entries:
        decomp,ehull = pd.get_decomp_and_e_above_hull(ent)
        form_atom = pd.get_form_energy_per_atom(ent)
        list_of_all.append([ent.composition.reduced_formula, float(form_atom), float(ehull)])
    

#    for i in list_of_all:
#        print(i)

    chosen_compounds = list()

    for i in list_of_all:
        if np.array(((i[2:]))).astype(np.float) <= e_above_hull:
#            print("Compound Matching criteria at ", i)
            chosen_compounds.append(i)
    
    return list_of_all, chosen_compounds
        
def extract_binaries(chosen_compounds):
    
    binaries_list = list()
    binaries_element_list = list()
    binaries_atomic_fraction = list()
    for i in chosen_compounds: 
        compound = mg.Composition(i[0])
        if len(compound.elements) == 2:
            binaries_list.append(i)
            binaries_element_list.append(compound.elements)
            fractional_comp = compound.fractional_composition
            
            for indx,i in enumerate(fractional_comp):
                if indx == 0:
                    x = compound.get_atomic_fraction(i)
                else:
                    y =  compound.get_atomic_fraction(i)
            
            binaries_atomic_fraction.append([x,y])
            
    return binaries_list, binaries_element_list, binaries_atomic_fraction

def build_binary_hull (binaries_list, binaries_element_list, binaries_atomic_fraction, element, index):


    chosen_binaries = list()
    chosen_elements = list()
    chosen_atomic_fractions = list()
    
    fraction_energy_list = list()
    
    for index, i in enumerate(index):
        if i == 1:
            chosen_binaries.append(binaries_list[index])
            chosen_atomic_fractions.append(binaries_atomic_fraction[index])
            chosen_elements.append(binaries_element_list[index])
    
    #print(chosen_binaries)
    #print(chosen_atomic_fractions)
    #print(chosen_elements)
    
    for index,i in enumerate(chosen_elements):
        if str(element) == str(i[0]):
            fraction_energy_list.append([chosen_atomic_fractions[index][0],chosen_binaries[index][1] ])
        elif str(element) == str(i[1]):
            fraction_energy_list.append([chosen_atomic_fractions[index][1],chosen_binaries[index][1] ])
        else:
            raise Exception('Element mis-matched. Check the inputs')
        
    
    return fraction_energy_list


def match_sides_from_binaries (binaries_element_list, three_elements):
    
    index_AC = [0]  * len(binaries_element_list)
    index_AB = [0]  * len(binaries_element_list)
    index_BC = [0]  * len(binaries_element_list)


    for index,i in enumerate(binaries_element_list):
        if str(three_elements[0]) == str(i[0]) and str(three_elements[2]) == str(i[1]) \
        or str(three_elements[0]) == str(i[1]) and str(three_elements[2]) == str(i[0]):
            index_AC[index] = 1   

        elif str(three_elements[0]) == str(i[0]) and str(three_elements[1]) == str(i[1]) \
        or str(three_elements[0]) == str(i[1]) and str(three_elements[1]) == str(i[0]):
            index_AB[index] = 1        

        elif str(three_elements[1]) == str(i[0]) and str(three_elements[2]) == str(i[1]) \
        or str(three_elements[1]) == str(i[1]) and str(three_elements[2]) == str(i[0]) :
            index_BC[index] = 1     
            
    return index_AC, index_AB, index_BC

def find_point(a,b,c, print_signal = 'ON'):
    xa = Symbol('xa')
    xc = Symbol('xc')

    xa1,xa2 = solve((1-(4*xa)+(4*(xa**2)))-a**2)
    if print_signal == 'ON': print('xa is', xa1,'or',xa2)

    xc1, xc2 = solve(((4*xc**2)-(8*xc)+4)-c**2)
    if print_signal == 'ON': print('xc is', xc1, 'or', xc2)

    ya1 = (math.sqrt(3)*xa1) 
    ya2 = (math.sqrt(3)*xa2) 
    if print_signal == 'ON': print('ya is', ya1, 'or', ya2)

    yc1 = -(math.sqrt(3)*xc1)+math.sqrt(3)
    yc2 = -(math.sqrt(3)*xc2)+math.sqrt(3)
    if print_signal == 'ON': print('yc is', yc1, 'or', yc2)


    m = Symbol('m')

    m1 = solve(((ya1 - yc1 + (math.sqrt(3)*xa1))/math.sqrt(3)) - m)
    n1 = yc1
    m2 = solve(((ya2 - yc2 + (math.sqrt(3)*xa2))/math.sqrt(3)) - m)
    n2 = yc2
    if print_signal == 'ON': print('m is', m1, 'or', m2)
    if print_signal == 'ON': print('n is', n1, 'or', n2)

    #choosing one of the two possible answers
    if print_signal == 'ON': print("Choosing one of the two possible answers...")
    if xa1 > 1 or xc1 > 1 or ya1 > 1 or yc1 > 1 or n1 > 1 or xa1 < 0 or xc1 < 0 or ya1 < 0 or yc1 < 0 or n1 < 0 :
        xa = xa2
        xc = xc2
        ya = ya2
        yc = yc2
        m = np.array(m2).astype(np.float)
        n = n2
    else :
        xa = xa1
        xc = xc1
        ya = ya1
        yc = yc1
        m = np.array(m1).astype(np.float)
        n = n1

    m = float(m)
    
    if print_signal == 'ON': 
        print('Final xa',xa)
        print('Final xc',xc)
        print('Final ya',ya)
        print('Final yc',yc)
        print('Final m',m)
        print('Final n',n)
    
    return [m,n]

def find_distance (x1,y1,x2,y2):
    return math.sqrt((x2-x1)**2+(y2-y1)**2)

def ratio_of_sides (m,n, print_signal = 'ON'):

    if print_signal == 'ON': print ('m,n',m,n)

    def equation2(p, *data):
        j, k = p
        m,n = data
        return ((j/math.sqrt(3))-(m/math.sqrt(3))-n+k, k-(j*math.sqrt(3)))

    data2 = (m,n)

    j, k =  fsolve(equation2,(1,1),args=data2)

    if print_signal == 'ON': print('j,k',j,k)
    

    
    aca = find_distance(0.5,math.sqrt(3)/2,j,k)
    acc = find_distance(0,0,j,k)
    total = aca + acc
    aca = aca / total
    acc = acc / total

    
    '''
    def equation3(p, *data):
        c1, c2 = p
        j,k = data
        return ((math.sqrt(j**2+k**2)/math.sqrt((0.5-j)**2+((math.sqrt(3)/2)-k)**2))-(c1/c2), c1+c2-1)

    data3 = (j,k)
    acc1, aca1 =  fsolve(equation3,(0.5,0.5),args=data3)
    
    
    if round(aca1,5) != round(aca,5):
        raise Exception('check formula for aca', aca, aca1)
    elif round(acc1,5) != round(acc,5):
        raise Exception('check formula for acc', acc, acc1)
    '''
    if print_signal == 'ON': 
        print('The triangle starts with A from the lower left corner and goes anticlockwise in the order of A-B-C.')
        print('acc,aca',acc,aca)


    abb = m
    aba = 1 - abb
    if print_signal == 'ON': print('abb,aba', abb, aba)

    def equation4(p, *data):
        j2, k2 = p
        m,n = data
        return (j2-m-(math.sqrt(3)*k2)+(math.sqrt(3)*n) , k2 + (math.sqrt(3)*j2) - math.sqrt(3))

    j2, k2 =  fsolve(equation4,(2,2),args=data2)

    if print_signal == 'ON': print('j2,k2',j2,k2)
   
    cbb = find_distance(0.5,math.sqrt(3)/2,j2,k2)
    cbc = find_distance(1,0,j2,k2)
    total = cbc + cbb
    cbc = cbc / total
    cbb = cbb / total
    
    if print_signal == 'ON': print('cbc,cbb',cbc,cbb)
    
    return [acc,aca,abb,aba,cbc,cbb]

def find_tri_composition_ratios (compound, print_signal = 'ON'):

    a = compound.get_atomic_fraction(compound.elements[0])
    b = compound.get_atomic_fraction(compound.elements[1])
    c = compound.get_atomic_fraction(compound.elements[2])
    
    if print_signal == 'ON': print('a,b,c,', a,b,c)
    return [a,b,c]

def CN_compounds (chargesof1st,chargesof2nd,chargesof3rd, max_nof_atom_per_compound=100, mix_charges='OFF', print_signal='ON'):
    
    combined_charges = chargesof1st + chargesof2nd + chargesof3rd
    max_nof_atom_per_element = max_nof_atom_per_compound-1
    total_num_of_charge_states = len(combined_charges)
        
    if mix_charges =='OFF':
        
        cn_atom_list = list()
        cn_charges_list = list()
        
        #mix-charge state - OFF
        for round1 in range(1,max_nof_atom_per_compound-1):
            for i in chargesof1st:
                for round2 in range(1,max_nof_atom_per_compound-1):
                    for j in chargesof2nd:
                        for round3 in range(1,max_nof_atom_per_compound-1):
                            for k in chargesof3rd:
                                addition = (i*round1) + (j*round2) + (k*round3)
                                if addition == 0 and (round1+round2+round3) < max_nof_atom_per_compound:
                                    cn_charges_list.append([i,j,k])
                                    cn_atom_list.append([round1,round2,round3])    
    
    
    else:
        

        no_atom_list = [ [0]* 3 for _ in range(pow(max_nof_atom_per_element,total_num_of_charge_states)) ]
        charges_array = [ [0]* (1+total_num_of_charge_states*(max_nof_atom_per_element-1)) for _ in range(pow(max_nof_atom_per_element,total_num_of_charge_states)) ]


        for indx1, i in enumerate(combined_charges):
            for r1 in range(pow(max_nof_atom_per_element,total_num_of_charge_states)):
                for fill_r1 in range(0,(int(r1/pow(max_nof_atom_per_element,indx1))%max_nof_atom_per_element)):
                    charges_array[r1][fill_r1+((max_nof_atom_per_element-1)*indx1)] = i



        for inde1, rows in enumerate(charges_array):

            count1 = 0;
            count2 = 0;
            count3 = 0;
            charges_sum = 0;

            for inde2,charges in enumerate(rows):
                charges_sum = charges_sum + charges
                
                #counting the number of zeros in the charges_array
                if inde2 < ((max_nof_atom_per_element-1)*len(chargesof1st)) and charges == 0:
                    count1 = count1 + 1
                elif inde2 < ((max_nof_atom_per_element-1)*(len(chargesof1st)+len(chargesof2nd))) and charges == 0:
                    count2 = count2 + 1
                elif inde2 < ((max_nof_atom_per_element-1)*(len(chargesof1st)+len(chargesof2nd)+len(chargesof3rd))) and charges == 0:
                    count3 = count3 + 1

            #the acttual number of atoms is the total minus the number of zeros
            no_atom_list[inde1][0] = ((max_nof_atom_per_element-1)*len(chargesof1st))-count1
            no_atom_list[inde1][1] = ((max_nof_atom_per_element-1)*len(chargesof2nd))-count2
            no_atom_list[inde1][2] = ((max_nof_atom_per_element-1)*len(chargesof3rd))-count3
            charges_array[inde1][inde2] = charges_sum    
            
        if print_signal == 'ON':
            print("Printing 5 rows")
            for indx, i in enumerate(charges_array):
                if indx < 5:
                    for j in i:
                        print(j, end = " ")
                    print()            
            #print(no_atom_list)

        cn_charges_list = list()
        cn_atom_list = list()

        
        for row_number,i in enumerate(charges_array):
            last_index = len(i)-1
            if i[last_index] == 0 and row_number > 0 and no_atom_list[row_number][0] > 0 and no_atom_list[row_number][1] > 0 and no_atom_list[row_number][2] > 0:
                #for j in i:
                #    print(j, end = " ")
                #print()
                cn_charges_list.append(i)
                cn_atom_list.append(no_atom_list[row_number])
        

    if print_signal == 'ON':
        print("Printing up to the first 10 charges")
        for ind, i in enumerate(cn_charges_list):
            if ind < 10:
                print(i)

        print("Printing the number of atoms")
        print(cn_atom_list)
    

    normalized_cn_atom_list = list()

    for ind, i in enumerate(cn_atom_list):
        n_row =([j/max(i) for j in i])
        divider = max(i)+1
        while sum(n_row)%1 != 0 or min(n_row) < 1:
            divider = divider - 1
            n_row =([j/divider for j in i])
        if sum([int(j/divider) for j in i]) == 0:
            normalized_cn_atom_list.append([int(j/j) for j in i])
        else:
            normalized_cn_atom_list.append([int(j/divider) for j in i])
    
    if print_signal == 'ON':
        print("Printing normalized atom list")
        print(normalized_cn_atom_list)
    
    return normalized_cn_atom_list, cn_charges_list       



def make_formula(three_elements,no_atoms):
    
    formula_list = list()
    element_list = list()
    
    for inx1, no_atoms in enumerate(no_atoms):
        formula_list.append(three_elements[0]+str(no_atoms[0])+three_elements[1]+str(no_atoms[1])+three_elements[2]+str(no_atoms[2]))
        element_list.append(three_elements)
        
    return formula_list, element_list

def get_energy_value (fraction_to_energy, ratio):
    
    sorted_list = [[0,0]]
    sorted_list.extend(sorted(fraction_to_energy, key=lambda x: x[0]))
    sorted_list.append([1,0])


    indx = 0
    while ratio > float(sorted_list[indx][0]):
        #print("Run")
        indx = indx + 1
        #print(indx)
        if indx == len(sorted_list)-1:
            break

    return ( ( sorted_list[indx][1] - sorted_list[indx-1][1] ) / ( sorted_list[indx][0] - sorted_list[indx-1][0] ) \
            * ( ratio - sorted_list[indx-1][0] ) ) +  sorted_list[indx-1][1] 


def get_breakdown_energy_values (energy1, ratio1, energy2, ratio2, energy3, ratio3, print_signal = 'ON'):

    r1 = get_energy_value(energy1, ratio1)
    r2 = get_energy_value(energy2, ratio2)
    r3 = get_energy_value(energy3, ratio3)
    
    if print_signal == 'ON': print('e1,e2,e3', r1,r2,r3)
    
    return [r1, r2, r3]

def find_weight_ratio(a,b,c, print_signal = 'ON'):
    
    weight = (1/a) + (1/b) + (1/c)
    
    #the corresponind values for w1 is not a becasue of inverse weight ratio
    w1 = (1/b)/weight
    w2 = (1/c)/weight
    w3 = (1/a)/weight
    
    if print_signal == 'ON': print('w1,w2,w3', w1,w2,w3)
    
    return w1,w2,w3

def find_lcm(a, b):
    return a * b // gcd(a, b)

def common_integer(numbers):
    
    #print('numbers', numbers)
    fractions = [Fraction(n).limit_denominator() for n in numbers]
    #print('fractions', fractions)
    multiple  = reduce(find_lcm, [f.denominator for f in fractions])
    #print('multiple', multiple)
    ints      = [f * multiple for f in fractions]
    #print('ints', ints)
    divisor   = reduce(gcd, ints)
    #print('divisor', divisor)
    return [int(n / divisor) for n in ints]


def get_atomic_composition_of_elements(ref, compound):

    original_comp = [0, 0, 0]
    
    ref = mg.Composition(ref)
    ref_elements = ref.elements
    #print(ref_elements)
    
    compound = mg.Composition(compound)
    factor = compound.get_reduced_formula_and_factor()[1]
    dict_val = compound.to_reduced_dict
    
    
    for indx, val in enumerate(ref_elements):
        original_comp[indx] = dict_val[str(val)]

    return [x*factor for x in original_comp]


def balance_the_equation(poriginal, pdecompo,iteration_limit=100):
    
    ref = poriginal
    original_compound = get_atomic_composition_of_elements(ref, poriginal)
    #print('original_compound',original_compound )
    
    
    decompo_compoounds = list()    
    for index, i in enumerate(pdecompo):
        compositions =[ x*(-1) for x in get_atomic_composition_of_elements(ref, pdecompo[index])]
        decompo_compoounds.append(compositions)
    #print('decompo_compoounds', decompo_compoounds)
    
    elementlist = ref
    elementMatrix = [original_compound]
    elementMatrix.extend(decompo_compoounds)
    
    elementMatrix = Matrix(elementMatrix)
    elementMatrix = elementMatrix.transpose()
    solution = elementMatrix.nullspace()[0]
    #multiplier = lcm([val for val in solution])
  
    nlist = list()
    for ind,i in enumerate([x for x in solution.tolist()]):
        nlist.extend(i)
        
    #print(elementlist)
    #print(elementMatrix)
    #print(nlist)
    multiplier = common_integer([str(x) for x in nlist])
    #print(multiplier)
    changed_format = [0, 0, 0, 0]
    for indx, val in enumerate(multiplier):
        changed_format[indx] = val+changed_format[indx]
    
    

    #print('changed_format', changed_format)
    
    
    
    '''
    result_right = [0, 0, 0]
    result_left = [0, 0, 1]
    multi = 1
    
    iteration = 0    
    
    

    while result_right != result_left:
        

        result_right = [0, 0, 0]
        

        for indx,val in enumerate(decompo_compoounds):
            result_tmp = [i*changed_format[indx+1]*(-multi)  for i in val]
            result_right = [a+b for a,b in zip(result_right,result_tmp)]

        result_left = [i*changed_format[0]*multi for i in original_compound]


        multi = multi + 1
        iteration = iteration + 1
        

        print('result_right', result_right)
        print('result_left', result_left)
        #print('m', multiplier)
        #print('l', result_left)
        #print('r', result_right)

        if iteration > iteration_limit:
            raise Exception('iteration_limit reached')
            break
    
    print('multi',multi-1)
    '''
    return tuple(changed_format)
    
    
def make_grid (three_elements, interval):
    no_atoms = list()
    format_no_atoms = list()
    for i in range(1,100,interval):
        for j in range(1,100,interval):
            if i+j < 100:
                no_atoms.append([i, j, 100-i-j])
                format_no_atoms.append((i,j,100-i-j))

    formula, element = make_formula(three_elements,no_atoms)
    
    return format_no_atoms, formula, element



def create_plotly_lines (z):
    
    minz = min(z)
    maxz = max(z)
    x=list([0, 1, 0.5, 0]*2)
    x.extend([0, 1, 1, 0.5, 0.5])
    y=list([0, 0, math.sqrt(3)/2, 0]*2)
    y.extend([0, 0,0,math.sqrt(3)/2,math.sqrt(3)/2])
    z = [minz,minz,minz,minz,maxz,maxz,maxz,maxz,minz,minz,maxz,maxz,minz]
    
    
    return go.Scatter3d( 
        x=x,
        y=y,
        z=z,
        mode="lines",
        hoverinfo="none",
        line=dict(color="rgba (0, 0, 0, 0.4)", dash="solid", width=1.0),
        showlegend=False,
    )

def create_plotly_mesh(m_t,n_t,z):
    
    x=[float(x) for x in m_t]
    y=[float(x) for x in n_t]
    z=[float(x) for x in z]
    
    return go.Mesh3d(
    x=x,
    y=y,
    z=z,
    opacity=0.5,
    colorscale = 'Viridis',
    colorbar=dict(title="Formation energy<br>(eV/atom)", x=0.9, len=0.75),
    intensity = z,
    hoverinfo="skip",
    lighting=dict(diffuse=0.0, ambient=1.0),
    flatshading=True,
    showlegend=True,
    #intensitymode='cell', 
)


def create_plotly_points(m_t,n_t,z,text, color, marker_size = 10):
    
    x=[float(x) for x in m_t]
    y=[float(x) for x in n_t]
    z=[float(x) for x in z]
    
    return go.Scatter3d(
    x=x,
    y=y,
    z=z,
    opacity=0.8,
    connectgaps = True,
    #hoverinfo="skip",
    mode="markers",
    marker = dict(color=color, size=marker_size, line=dict(color="black", width=2)),
    hoverinfo = "all", 
    hovertext = text,
    #marker=go.scatter3d.Marker(color=z, size=marker_size, showscale=False, colorscale='Viridis', colorbar=dict(thickness=10))
    )


def get_default_annotation ():
    return dict(
            showarrow=False,
            x=0,
            y=0,
            z=0,
            text="Point 1",
            xanchor="left",
            xshift=10,
            font = {'color': 'black', 'size': 25},
            opacity=0.9,) 

def get_annotations_list(three_elements):
    annotations_list = list()

    for i in [0,1,2]:
        annotation = get_default_annotation()
        if i == 0:
            annotation.update({"x" : 0,
                              "y" : 0,
                              "z" : -0.01,
                              "text" : three_elements[0]})
        elif i == 1:
            annotation.update({"x" : 1,
                              "y" : 0,
                              "z" : -0.01,
                              "text" : three_elements[1]})
        elif i == 2:
            annotation.update({"x" : 0.5,
                              "y" : math.sqrt(3)/2,
                              "z" : -0.01,
                              "text" : three_elements[2],
                              "xanchor" : "right",
                              "xshift" : -10})
        annotations_list.append(annotation)
    return annotations_list

def get_default_layout (three_elements):
    layout = {}
    layout.update({'autosize': True,
    'height': 700,
    'hovermode': 'closest',
    'paper_bgcolor': 'rgba(0,0,0,0)',
    'plot_bgcolor': 'rgba(0,0,0,0)',
    'margin': {'b': 10, 'l': 0, 'pad': 0, 't': 0, 'r': 0},
    'showlegend': True,
    'legend': {'orientation': 'h',
    'x': 0.5,
    'y': 0.0,
    'traceorder': 'reversed',
    'xanchor': 'center',
    'yanchor': 'top'},
    'scene_camera': {'center': {'x': -0.1, 'y': 0, 'z': -0.15},
    'eye': {'x': -0.1, 'y': 0, 'z': 2.5},
    'projection': {'type': 'orthographic'}},
    'scene': {'xaxis': {'title': None,
    'visible': False,
    'autorange': True,
    'showgrid': False,
    'zeroline': False,
    'showline': False,
    'ticks': '',
    'showaxeslabels': False,
    'showticklabels': False,
    'showspikes': False},
    'yaxis': {'title': None,
    'visible': False,
    'autorange': True,
    'showgrid': False,
    'zeroline': False,
    'showline': False,
    'ticks': '',
    'showaxeslabels': False,
    'showticklabels': False,
    'showspikes': False},
    'zaxis': {'title': None,
    'visible': False,
    'autorange': True,
    'showgrid': False,
    'zeroline': False,
    'showline': False,
    'ticks': '',
    'showaxeslabels': False,
    'showticklabels': False,
    'showspikes': False},
    'annotations': [{'showarrow': False,
    'x': 0,
    'y': 0,
    'z': -0.01,
    'text': 'Cu',
    'xanchor': 'left',
    'xshift': 10,
    'font': {'color': 'black', 'size': 25},
    'opacity': 0.9},
    {'showarrow': False,
    'x': 1,
    'y': 0,
    'z': -0.01,
    'text': 'Bi',
    'xanchor': 'left',
    'xshift': 10,
    'font': {'color': 'black', 'size': 25},
    'opacity': 0.9},
    {'showarrow': False,
    'x': 0.5,
    'y': 0.8660254037844386,
    'z': -0.01,
    'text': 'S',
    'xanchor': 'right',
    'xshift': -10,
    'font': {'color': 'black', 'size': 25},
    'opacity': 0.9}]},
    'scene_aspectratio': {'x': 1.7, 'y': 1.7, 'z': 1.2}})
    
    annotations_list = get_annotations_list(three_elements)
    layout["scene"].update({"annotations": annotations_list})
    
    return layout

def generate_ploty_data(m_t,n_t,z_t,CN_pretty_formula,three_elements, marker_size=10):
    
    data = create_plotly_mesh(m_t,n_t,z_t)
    lines = create_plotly_lines(z_t)
    points = create_plotly_points(m_t[-len(CN_pretty_formula):],n_t[-len(CN_pretty_formula):],z_t[-len(CN_pretty_formula):],CN_pretty_formula,"darkgreen", marker_size)
    
    return [data, lines, points]




def generate_ternary_plot (no_atoms_inc, data, title, three_elements, scale = 100):

    pe_dic = {}
    for i in range(len(no_atoms_inc)):
        pe_dic[no_atoms_inc[i]] = data[i]

    figure, tax = ternary.figure(scale=scale)
    figure.set_size_inches(25, 16)
    tax.heatmap(pe_dic, scale=scale, style="h")
    tax.boundary(linewidth=3.0)
    tax.set_title(three_elements[0]+' '+three_elements[1]+' '+three_elements[2]+' '+str(title) , fontsize=25)
    tax.ticks(axis='lbr', linewidth=1, multiple=10)
    tax.left_axis_label(three_elements[2],fontsize=30)
    tax.right_axis_label(three_elements[1],fontsize=30)
    tax.bottom_axis_label(three_elements[0],fontsize=30)
    tax.clear_matplotlib_ticks()
    
    return figure, tax

def triangular_coord(coord):
    """
    Convert a 2D coordinate into a triangle-based coordinate system for a
    prettier phase diagram.
    Args:
        coord: coordinate used in the convex hull computation.
    Returns:
        coordinates in a triangular-based coordinate system.
    """
    unitvec = np.array([[1, 0], [0.5, math.sqrt(3) / 2]])

    result = np.dot(np.array(coord), unitvec)
    return result.transpose()

def area(x1, y1, x2, y2, x3, y3):
 
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0)

def isInside(x1, y1, x2, y2, x3, y3, x, y):
 
    # Calculate area of triangle ABC
    A = area (x1, y1, x2, y2, x3, y3)
 
    # Calculate area of triangle PBC, PAC, PAB
    A1 = area (x, y, x2, y2, x3, y3)
    A2 = area (x1, y1, x, y, x3, y3)
    A3 = area (x1, y1, x2, y2, x, y)
     
    # Check if sum of A1, A2 and A3 is same as A
    #print(x1, y1, x2, y2, x3, y3, x, y)
    #print('A', A)
    #print('A1+A2+A3', A1, A2, A3, A1+A2+A3)
    if(A + 1E-15  > A1 + A2 + A3):
        #print('True')
        return True
    else:
        #print('False')
        return False
    
def define_area(point1, point2, point3):

    point1 = np.asarray(point1)
    point2 = np.asarray(point2)
    point3 = np.asarray(point3)
    AB = np.asmatrix(point2 - point1)
    AC = np.asmatrix(point3 - point1)
    N = np.cross(AB, AC) # Vector cross product, find the normal vector
    
    # Ax+By+Cz+D = 0
    Ax = N[0, 0]
    By = N[0, 1]
    Cz = N[0, 2]
    D = -(Ax * point1[0] + By * point1[1] + Cz * point1[2])
    return Ax, By, Cz, D

def point2area_distance(Ax, By, Cz, D, point4):

    #mod_d = Ax * point4[0] + By * point4[1] + Cz * point4[2] + D
    #mod_area = np.sqrt(np.sum(np.square([Ax, By, Cz])))
    #d = abs(mod_d) / mod_area
    
    pointZ = ( -D - Ax * point4[0] - By * point4[1]) / Cz
    
    d = pointZ - point4[2]

    return pointZ, d

def get_distance_from_hull(list_of_vertices,m,n, z_energy, print_signal = 'ON'):
    
    if print_signal == 'ON': print("\nFinding distance from MP's hull surface...")
       
    triangle_indicator = get_triangle_indicator(list_of_vertices,m,n)
    
    base_triangle = get_base_triangle (list_of_vertices, triangle_indicator, print_signal)
    
        
    point4 = [m, n, z_energy]
    Ax, By, Cz, D = base_triangle
    pointZ, d1 = point2area_distance(Ax, By, Cz, D, point4)
    
    if print_signal == 'ON': print("Distance from hull is ", d1)
    
    return pointZ, d1, triangle_indicator
    
def get_triangle_indicator(list_of_vertices,m,n):

    triangle_indicator = -1
    
    for indx, valx in enumerate(list_of_vertices):
        x1 = valx[0][0]
        y1 = valx[0][1]
        x2 = valx[1][0]
        y2 = valx[1][1]
        x3 = valx[2][0]
        y3 = valx[2][1]

        if isInside(x1, y1, x2, y2, x3, y3, m, n):
            triangle_indicator = indx

    if triangle_indicator == -1:
        raise Exception('Point (m,n) is not within the given list_of_vertices')
    else:
        return triangle_indicator
    
def get_base_triangle (list_of_vertices, triangle_indicator, print_signal):
    
    point1 = list_of_vertices[triangle_indicator][0]
    point2 = list_of_vertices[triangle_indicator][1]
    point3 = list_of_vertices[triangle_indicator][2]
    Ax, By, Cz, D = define_area(point1, point2, point3)
        
    if print_signal == 'ON': print("Point is in the triangle with vertices of ", point1[0:2], point2[0:2], point3[0:2])
        
    return Ax, By, Cz, D

def get_list_of_vertices(pd):

    facets = np.array(pd.facets)
    coords = np.array([triangular_coord(c) for c in zip(pd.qhull_data[:-1, 0], pd.qhull_data[:-1, 1])])
    energies = np.array([pd.get_form_energy_per_atom(e) for e in pd.qhull_entries])
    
   
    list_of_vertices = list()
    for indx, valx in enumerate(facets):
        vertices = list()
        for val in (facets[indx]):
            #print(val)
            vertice = [coords[:, 0][val], coords[:, 1][val], energies[val] ]
            #print(vertice)
            vertices.append(vertice)

        list_of_vertices.append(vertices)
    
    return list_of_vertices

def find_decomposition_details(compound_name, key, pd = 0, list_of_vertices = 0, print_signal = 'ON'):
    
    compound = mg.Composition(compound_name)
    three_elements_fc = [str(compound.elements[0]), str(compound.elements[1]), str(compound.elements[2])]
    
    if pd == 0:
        if print_signal == 'ON': print('Generating new pd')
        pd = get_phasediagram(key, three_elements_fc)
        list_of_vertices = get_list_of_vertices(pd)
    
    #finding the decomposition products if it is just 1 decomposition subset, skip this compound
    decomposition_details = pd.get_decomposition(compound)
    decompos = [k.composition.reduced_formula for k in decomposition_details]
    
    if list_of_vertices == 0:
        raise Exception('Check the input list_of_vertices and most likely pd as well')
        
    return compound, decompos, pd, three_elements_fc, list_of_vertices

def analyze_compound (pd, three_elements, decompos, compound,e_above_hull, print_signal = 'ON'):
    compund_energy = [0, 0, 0, 0]

    [list_of_all, chosen_compounds] = get_specific_compounds (pd, e_above_hull)
    binaries_list, binaries_element_list, binaries_atomic_fraction = extract_binaries(chosen_compounds)
    index_AC, index_AB, index_BC = match_sides_from_binaries (binaries_element_list, three_elements)
    fraction_to_energy1 = build_binary_hull(binaries_list,binaries_element_list, binaries_atomic_fraction, three_elements[0], index_AC)
    fraction_to_energy2 = build_binary_hull(binaries_list,binaries_element_list, binaries_atomic_fraction, three_elements[1], index_AB)
    fraction_to_energy3 = build_binary_hull(binaries_list,binaries_element_list, binaries_atomic_fraction, three_elements[2], index_BC)

    #find energy of the decomposition compounds
    for indx, valx in enumerate(decompos):
        for indy, valy in enumerate(chosen_compounds):
            if str(valx) == str(valy[0]):
                de_comp = mg.Composition(valx)
                compund_energy[indx+1] = de_comp.num_atoms*valy[1]

    #find composiiton
    if print_signal == 'ON': print("\nFinding Composition...")
    a,b,c = find_tri_composition_ratios(compound, print_signal)

    #find the point in phasediagram
    if print_signal == 'ON': print("\nFinding point in phasediagram...")
    m,n = triangular_coord([a,c])

    #find the ratio of sides
    if print_signal == 'ON': print("\nFinding the ratio of sides...")
    acc,aca,abb,aba,cbc,cbb = ratio_of_sides(m,n, print_signal)

    #find the energy of the corresponding binaries
    if print_signal == 'ON': print("\nFinding the corresponding energy on each sides")
    e1, e2, e3 = get_breakdown_energy_values(fraction_to_energy1, aca, fraction_to_energy2, abb, fraction_to_energy3, cbc, print_signal)
   
    
    #check if the order of element is consistent    
    if str(three_elements[0]) != str(compound.elements[0])  \
    and str(three_elements[1]) != str(compound.elements[1]) \
    and str(three_elements[2]) != str(compound.elements[2]) :
        print(three_elements, compund.elements)
        raise Exception('Element mis-matched. Check the inputs')

    #find the weight ratio of energies
    if print_signal == 'ON': print("\nFinding the weight ratio of energies")
    w1, w2, w3 = find_weight_ratio(a,b,c, print_signal)

    #predict the formation energy per atom of the compund
    if print_signal == 'ON': print("\nFinding the formation energy of the compounds")
    compound_energy_per_atom =  (w1*e1) + (w2*e2) + (w3*e3)
    compund_energy[0] = compound_energy_per_atom * compound.num_atoms
    if print_signal == 'ON': print('compound_energy_per_atom', compound_energy_per_atom)
    if print_signal == 'ON': print('compound_energy', compund_energy)


    #balance the breakdown equation
    if print_signal == 'ON': print("\nBalancing the breakdown equation")
    multiplier = balance_the_equation(compound, decompos)
    if print_signal == 'ON': print(multiplier)


    #find the energy difference **note positive result means original has higher energy
    if print_signal == 'ON': print("\nFinding the energy difference between original and breakdown subsets")
    result = [i*y for i,y in zip(multiplier,compund_energy)]
    result_left = result[0]
    result_right = result[1] + result[2] + result[3]
    difference = result_left - result_right
    difference_per_atom = difference/compound.num_atoms/multiplier[0]
    if print_signal == 'ON': print('Left side ', result_left)
    if print_signal == 'ON': print('Right side ', result_right)
    if print_signal == 'ON': print('Energy difference is ', difference, ' eV.')
    if print_signal == 'ON': print('Energy difference/atom is ', difference_per_atom, ' eV/atom.')

    return m,n, compound_energy_per_atom, compund_energy, difference_per_atom


def get_hull_mesh (pd, colorscale, colorbar_title, shading_title, intensitymode):

    facets = np.array(pd.facets)
    coords = np.array([triangular_coord(c) for c in zip(pd.qhull_data[:-1, 0], pd.qhull_data[:-1, 1])])
    energies = np.array([pd.get_form_energy_per_atom(e) for e in pd.qhull_entries])

    x=list(coords[:, 0])
    y=list(coords[:, 1])
    z=list(energies)
    i=list(facets[:, 1])
    j=list(facets[:, 0])
    k=list(facets[:, 2])
    
    return go.Mesh3d(
            x=x,
            y=y,
            z=z,
            i=i,
            j=j,
            k=k,
            opacity=0.8,
            intensity=list(energies),
            colorscale=colorscale,
            colorbar=dict(title=colorbar_title, x=0.9, len=0.75),
            hoverinfo="none",
            lighting=dict(diffuse=0.0, ambient=1.0),
            name=shading_title,
            flatshading=True,
            showlegend=True,
            intensitymode=intensitymode #'cell' or 'vertex' (default)
        )

def main (formula_list, three_elements_list, key, e_above_hull, print_signal):
    #initialization of parameters
    num_decomposition_subset_list = list() #the number of decomposition products 
    m_t = list() #coordinates for point in phase diagram 
    n_t = list() #coordinate for point in phase diagram
    compound_t = list() #name of the compound
    compound_energy_per_atom_t = list() #formation energy of the compound per atom (eV/atom)
    compound_energy_t = list() #formation energy of the compound (eV)
    energy_difference_per_atom_t = list() #difference between original and breakdown product (FEPA)
    distance_from_hull_t = list() #the energy above MP hull to the calculated FE
    indexing = [0]*len(formula_list)
    triangle_indicator = list()

    for ind, val in enumerate(formula_list):
        compound_name = val
        three_elements = three_elements_list[ind]
    
        #if two compunds have the same phase diagram, this code enquire the pd only once to save time. compound list should be in alphabetical order
        if ind == 0:
            compound, decompos, pd, three_elements_fc,list_of_vertices = find_decomposition_details(compound_name, key, print_signal = print_signal)
        elif three_elements == three_elements_list[ind-1]:
            compound, decompos, pd, three_elements_fc,list_of_vertices = find_decomposition_details(compound_name, key,pd,list_of_vertices, print_signal = print_signal)
        else:
            compound, decompos, pd, three_elements_fc,list_of_vertices = find_decomposition_details(compound_name, key, print_signal = print_signal)

        if print_signal == 'ON': print(compound, ' will decompose into ', len(decompos), ' subsets:',decompos)
        num_decomposition_subset_list.append(decompos)

        #if there is only one decomposition product, there was no decomposition
        if len(decompos) > 1:

            #calculating point in phase diagram, FE of the compund it's decompositon products
            m,n, compound_energy_per_atom, compund_energy, difference_per_atom = analyze_compound (pd,three_elements_fc, decompos, compound, e_above_hull, print_signal)
    
            #getting the hull from MP and calculating the distance from hull    
            distance_from_hull, t_indicator = get_distance_from_hull(list_of_vertices, m, n, compound_energy_per_atom, print_signal)
        
        
            #append the results to the list of variables
            m_t.append(m)
            n_t.append(n)
            compound_t.append(compound_name)
            compound_energy_per_atom_t.append(compound_energy_per_atom)
            compound_energy_t.append(compund_energy)
            energy_difference_per_atom_t.append(difference_per_atom)
            distance_from_hull_t.append(distance_from_hull)
            indexing[ind] = 1
            triangle_indicator.append(t_indicator)

            if print_signal == 'ON': print('\n\n-------------------\n\n-------------------')

    return m_t,n_t,compound_t,compound_energy_per_atom_t,compound_energy_t,energy_difference_per_atom_t,distance_from_hull_t,indexing,triangle_indicator, pd

def remove_compound_duplicates (pretty_formula_list, three_elements_list,icsd_id, df_index, print_signal = 'ON'):
    
    n_icsd_id_list = list()
    n_pretty_formula_list = list()
    n_three_elements_list = list()
    n_df_index = list()

    #the first filter is removing those duplicates with no ICSD_ID, if there are two items with ICSD_ID, the first one is recorded.
    for ind, val in enumerate (pretty_formula_list):
        #print(ind)
        if ind > 0:
            if val == pretty_formula_list[ind-1] and len(n_icsd_id_list[-1]) != 0:
                if print_signal == 'ON' : print('I1')

            elif val == pretty_formula_list[ind-1] and len(n_icsd_id_list[-1]) == 0:

                if print_signal == 'ON' : print('I3')
                n_icsd_id_list.pop()
                n_pretty_formula_list.pop()
                n_three_elements_list.pop()
                n_df_index.pop()

                n_icsd_id_list.append(icsd_id[ind])
                n_pretty_formula_list.append(pretty_formula_list[ind])
                n_three_elements_list.append(three_elements_list[ind])
                n_df_index.append(df_index[ind])


            else:
                if print_signal == 'ON' : print('I2')
                n_icsd_id_list.append(icsd_id[ind])
                n_pretty_formula_list.append(pretty_formula_list[ind])
                n_three_elements_list.append(three_elements_list[ind])
                n_df_index.append(df_index[ind])

        else:
            n_icsd_id_list.append(icsd_id[ind])
            n_pretty_formula_list.append(pretty_formula_list[ind])
            n_three_elements_list.append(three_elements_list[ind])
            n_df_index.append(df_index[ind])

    return n_pretty_formula_list, n_three_elements_list, n_icsd_id_list, n_df_index
