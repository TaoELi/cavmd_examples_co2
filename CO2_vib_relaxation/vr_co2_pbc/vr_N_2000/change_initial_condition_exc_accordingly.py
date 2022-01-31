import numpy as np
import sys
import xml.etree.ElementTree as ET
import json as js
from collections.abc import Mapping, Sequence

NUM_ATOM_EACH_MOLECULE = 3
QUANTUM = 2000 # in units of K
FINAL_ENERGY = QUANTUM * 3 # in units of K
FINAL_ENERGY_FLUCT = QUANTUM * 0.25
N_MODIFIED = int(10.0 / 216.0 * 2000.0)
print("We will excite %d molecules" %N_MODIFIED)

class ScientificNotationEncoder(js.JSONEncoder):
    def iterencode(self, o, _one_shot=False):
        if isinstance(o, float):
            return "{:e}".format(o)
        elif isinstance(o, Mapping):
            return "{{{}}}".format(', '.join('"{}" : {}'.format(str(ok), self.iterencode(ov))
                                             for ok, ov in o.items()))
        elif isinstance(o, Sequence) and not isinstance(o, str):
            return "[{}]".format(', '.join(map(self.iterencode, o)))
        return ', '.join(super().iterencode(o, _one_shot))

def increase_kinetic_energy(p, m, which_molecule=0, num_atom=NUM_ATOM_EACH_MOLECULE, final_energy = FINAL_ENERGY):
    # Each molecule contains three atoms, here we increase the kinetic energy of one molecule
    # First, we obtain the index of the molecule we are interested in
    p_index = [which_molecule*num_atom*3 + i for i in range(num_atom*3)]
    #print("Momenta index for No. %d molecule is" %which_molecule, p_index)
    m_index = [which_molecule*num_atom + i for i in range(num_atom)]
    #print("Mass index for No. %d molecule is" %which_molecule, m_index)
    # Second, we obtain the momenta and mass for this molecule
    p_array = np.array(p)
    m_arraay = np.array(m)
    p_onemolecule = p_array[p_index]
    m_onemolecule = m_arraay[m_index]
    #print("Momenta for this molecule is", p_onemolecule)
    #print("Mass for this molecule is", m_onemolecule)
    # Third, we calculate the original kinetic energy
    #print(np.reshape(p_onemolecule, (-1, num_atom)))
    data = np.reshape(p_onemolecule**2, (-1, num_atom)) / m_onemolecule.reshape((-1, 1)) / 2.0
    DoFs = 3.0 * num_atom
    KE = np.sum(data) / DoFs * 315777.09
    #print("previous kinetic energy of molecule %d is %.2f K" %(which_molecule, KE))
    # Forth, we increas the kinetic energy by multiplying a factor to momenta
    # The final energy needs to be resampled with a random fluctuation

    # Method I:
    #final_energy += FINAL_ENERGY_FLUCT * (np.random.rand()-0.5)
    #factor = (final_energy / KE)**(0.5)
    #p_onemolecule = p_onemolecule * factor

    # method II:
    data_modified = (FINAL_ENERGY_FLUCT * (np.random.rand(3, 3)- 0.5) + final_energy) / 315777.09
    p_onemolecule_new = np.reshape(np.sqrt(data_modified * m_onemolecule.reshape((-1, 1)) * 2.0), -1)
    # random generation
    tominus_idx = np.random.rand(np.size(data_modified)) > 0.5
    p_onemolecule_new[tominus_idx] *= -1.0
    p_onemolecule = p_onemolecule_new

    data = np.reshape(p_onemolecule**2, (-1, num_atom)) / m_onemolecule.reshape((-1, 1)) / 2.0
    KE_new = np.sum(data) / DoFs * 315777.09
    #print("new kinetic energy of molecule %d is %.2f K" %(which_molecule, KE_new))
    p_array[p_index] = p_onemolecule
    return p_array.tolist(), KE, KE_new


def operate_one_file(xml_path, store_folder="modified_configs_exc_accordingly/"):
    # Open original file
    tree = ET.parse(xml_path)
    # Read this file
    root = tree.getroot()
    a = root.find("system")
    suba = a.find("beads")
    # load mass and momenta
    momenta = suba.find("p")
    mass = suba.find("m")
    momenta_value = momenta.text
    mass_value = mass.text
    momenta_value = js.loads(momenta_value)
    mass_value = js.loads(mass_value)
    # We define which molecules need to be modified from a list
    molecule_indices = np.arange(N_MODIFIED)
    # change the momenta and increase the kinetic energy
    ke_old_lst, ke_new_lst = [], []
    for i in molecule_indices:
        momenta_value, ke_old, ke_new = increase_kinetic_energy(momenta_value, mass_value,
                    which_molecule=i, final_energy = FINAL_ENERGY)
        ke_old_lst.append(ke_old)
        ke_new_lst.append(ke_new)

    ke_old_lst = np.array(ke_old_lst)
    ke_new_lst = np.array(ke_new_lst)

    print("Statistics: previous avg KE = %.2f K, std = %.2f K" %(np.mean(ke_old_lst), np.std(ke_old_lst)) )
    print("Statistics: new avg KE = %.2f K, std = %.2f K" %(np.mean(ke_new_lst), np.std(ke_new_lst)) )

    momenta.text = js.dumps(momenta_value, cls=ScientificNotationEncoder)

    # Write back to file
    pure_path = xml_path.split('/')[-1].split('_')[-1]
    tree.write(store_folder+"/init_" + pure_path)

    # Also write useful information about the reset of initial condition
    data = np.zeros((N_MODIFIED, 3))
    data[:, 0] = molecule_indices
    data[:, 1] = ke_old_lst
    data[:, 2] = ke_new_lst
    data_filename = store_folder+"/init_"+pure_path.split(".")[0]+".md"
    np.savetxt(data_filename, data, header="indices, old avergaed KE [K], new averaged KE [K]", fmt='%.4E')



if __name__ == "__main__":
    for path in sys.argv[1:]:
        print("Dealing with %s" %path)
        operate_one_file(xml_path=path)
