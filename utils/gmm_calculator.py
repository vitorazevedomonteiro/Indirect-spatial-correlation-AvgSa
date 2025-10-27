import numpy as np
import pandas as pd
from pathlib import Path
from openquake.hazardlib.gsim.aristeidou_2024 import AristeidouEtAl2024
from openquake.hazardlib.contexts import full_context, simple_cmaker
from openquake.hazardlib.contexts import SitesContext, RuptureContext, DistancesContext
from openquake.hazardlib import imt

def compute_gmm_predictions(database_path: Path, output_dir: Path, avgsa_periods):
    data = pd.read_csv(database_path)
    gmm = AristeidouEtAl2024()
    output_dir.mkdir(parents=True, exist_ok=True)

    vs30 = np.array(data['Vs30'])
    Mw = np.array(data['magnitude'])
    rake = np.array(data['rake'])
    Z2pt5 = np.array(data['Z2pt5'])
    Rrup = np.array(data['Rrup'])
    ztor = np.array(data['Ztor'])
    hypo_depth = np.array(data['D_hyp'])
    rjb = np.array(data['Rjb'])
    rx = np.array(data['Rx'])
    RSNs = data['RSN'].values
    EQIDs = data['EQID'].values
    sids = np.arange(len(vs30))

    for T in avgsa_periods:
        periods = np.linspace(0.2 * T, 2.0 * T, 10)
        im_objs = [imt.SA(p) for p in periods]
        results = []

        for idx in range(len(data)):
            sites = SitesContext()
            sites.__dict__.update({'z2pt5': Z2pt5[idx], 'vs30': vs30[idx], 'sids': np.array([sids[idx]])})
            rup = RuptureContext()
            rup.__dict__.update({'rake': rake[idx], 'mag': Mw[idx], 'ztor': ztor[idx], 'hypo_depth': hypo_depth[idx]})
            dists = DistancesContext()
            dists.__dict__.update({'rrup': Rrup[idx], 'rjb': rjb[idx], 'rx': rx[idx]})

            ctx = full_context(sites, rup, dists)
            cmaker = simple_cmaker([gmm], [obj.string for obj in im_objs], mags=[str(Mw[idx])])
            ctx = cmaker.recarray([ctx])

            mean = np.zeros((len(im_objs), 1))
            sig = np.zeros((len(im_objs), 1))
            tau = np.zeros((len(im_objs), 1))
            phi = np.zeros((len(im_objs), 1))
            gmm.compute(ctx, im_objs, mean, sig, tau, phi)

            for i, im_obj in enumerate(im_objs):
                results.append({
                    'RSN': RSNs[idx],
                    'EQID': EQIDs[idx],
                    'Period': im_obj.period,
                    'Mean': np.exp(mean[i][0]),
                    'Stdev1': sig[i][0],
                    'Stdev2': tau[i][0],
                    'Stdev3': phi[i][0]
                })

        output_file = output_dir / f"predicted_AvgSa2_sa({T:.2f}).csv"
        pd.DataFrame(results).to_csv(output_file, index=False)
        print(f"Saved: {output_file}")
