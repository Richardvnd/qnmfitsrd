import pickle 
import json
import qnmfitsrd as qnmfits

filepath  = '/data/rvnd2-2/CCE_data/superrest_data'

def SXS_CCE(ID, zero_time=(2,2), lev = 'Lev5', radius = 'R2'):
        
    if ID == '0305':

        print("Note that 0305 only has one level and radius. These parameters will be ignored.")
        
        with open(f'{filepath}/SXS:BBH_ExtCCE_superrest:{ID}/SXS:BBH_ExtCCE_superrest:{ID}.pickle', 'rb') as f:
            h_prime_dict = pickle.load(f)
        with open(f'{filepath}/SXS:BBH_ExtCCE_superrest:{ID}/SXS:BBH_ExtCCE_superrest:{ID}_metadata.json', 'r') as f:
            metadata = json.load(f)
    
    else:
    
        with open(f'{filepath}/SXS:BBH_ExtCCE_superrest:{ID}/SXS:BBH_ExtCCE_superrest:{ID}_{lev}_{radius}.pickle', 'rb') as f:
            h_prime_dict = pickle.load(f)
        with open(f'{filepath}/SXS:BBH_ExtCCE_superrest:{ID}/SXS:BBH_ExtCCE_superrest:{ID}_{lev}_{radius}_metadata.json', 'r') as f:
            metadata = json.load(f)

    times = h_prime_dict.pop('times')
    
    sim = qnmfits.Custom(
            times, 
            h_prime_dict, 
            metadata={'remnant_mass':metadata['M_f'], 'remnant_dimensionless_spin':metadata['chi_f']}, 
            zero_time=zero_time
            )
    
    return sim 