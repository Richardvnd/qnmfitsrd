import sxs  
import scri  
import spherical_functions as sf  
import numpy as np
import pickle
import matplotlib.pyplot as plt
import json 

from scri.asymptotic_bondi_data.map_to_superrest_frame import MT_to_WM

# List of second smallest radii
filenum = ['0002', '0003', '0004', '0005', '0006', '0007', '0008', '0009', '0010', '0011', '0012', '0013']
radii = ['0261', '0250', '0236', '0274', '0273', '0270', '0305', '0270', '0235', '0222', '0223', '0237']

for filenum, radius in zip(filenum, radii):

    filepath = f"/data/rvnd2-2/CCE_data"
    filepathraw = f"{filepath}/raw_data/SXS:BBH_ExtCCE:{filenum}"

    abd = scri.SpEC.file_io.create_abd_from_h5(
        h=f"{filepathraw}/Lev5:rhOverM_BondiCce_R{radius}.h5",
        Psi4=f"{filepathraw}/Lev5:rMPsi4_BondiCce_R{radius}.h5",
        Psi3=f"{filepathraw}/Lev5:r2Psi3_BondiCce_R{radius}.h5",
        Psi2=f"{filepathraw}/Lev5:r3Psi2OverM_BondiCce_R{radius}.h5",
        Psi1=f"{filepathraw}/Lev5:r4Psi1OverM2_BondiCce_R{radius}.h5",
        Psi0=f"{filepathraw}/Lev5:r5Psi0OverM3_BondiCce_R{radius}.h5",
        file_format="RPDMB",
    )

    # Define t = 0 to be the time of peak luminosity
    abd.t -= abd.t[np.argmax(MT_to_WM(2.0 * abd.sigma.bar.dot).norm())]

    # Interpolate to the merger/ringdown regime to speed up the fixing of the BMS frame
    abd_ringdown = abd.interpolate(np.arange(-100, abd.t[-1], 0.1))

    # Map to the remnant superrest frame (pick some time window after QNMs have decayed - 300 recommended by Eliot)
    abd_ringdown_superrest, _, _ = abd_ringdown.map_to_superrest_frame(t_0=300, padding_time=20)

    # Compute the remnant mass and spin
    M_f = abd_ringdown_superrest.bondi_rest_mass()[-1]
    chi_f = abd_ringdown_superrest.bondi_dimensionless_spin()[-1]
    chi_f = chi_f.tolist() 

    metadata = {'M_f': M_f, 'chi_f': chi_f}

    # Compute the strain waveform
    # You can save/load this with, e.g., scri.SpEC.file_io.write_to_h5/read_from_h5
    h = MT_to_WM(2.0 * abd_ringdown_superrest.sigma.bar)

    h_dict = {'times': h.t}

    for ell in range(2, h.ell_max+1):
        for m in range(-ell,ell+1):
            h_dict[ell,m] = np.array(h.data[:, h.index(ell,m)])

    with open(f'{filepath}/superrest_data/SXS:BBH_ExctCCE:{filenum}_superrest/SXS:BBH_ExctCCE:{filenum}_superrest.pkl', 'wb') as f:
        pickle.dump(h_dict, f)

    with open(f'{filepath}/superrest_data/SXS:BBH_ExctCCE:{filenum}_superrest/SXS:BBH_ExctCCE:{filenum}_superrest_metadata.json', 'w') as f:
        json.dump(metadata, f)