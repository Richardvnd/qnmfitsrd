import pickle 
import json
import qnmfitsrd as qnmfits

filepath  = '/data/rvnd2-2/CCE_data/superrest_data'

def SXS_CCE(ID, zero_time=(2,2)):
    with open(f'{filepath}/SXS:BBH_ExctCCE:{ID}_superrest/SXS:BBH_ExctCCE:{ID}_superrest.pkl', 'rb') as f:
        h_prime_dict = pickle.load(f)
    with open(f'{filepath}/SXS:BBH_ExctCCE:{ID}_superrest/SXS:BBH_ExctCCE:{ID}_superrest_metadata.json', 'r') as f:
        metadata = json.load(f)

    times = h_prime_dict.pop('times')
    
    sim = qnmfits.Custom(
    times, 
    h_prime_dict, 
    metadata={'remnant_mass':metadata['M_f'], 'remnant_dimensionless_spin':metadata['chi_f']}, 
    zero_time=zero_time
    )

    return sim