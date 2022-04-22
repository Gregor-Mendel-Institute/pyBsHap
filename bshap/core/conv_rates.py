
import pyro
import pyro.distributions as dist
import torch
import numpy as np


class CytosineModel:
    ## has three parameters
    # 1. Depth
    #    number of reads you sample
    # 2. relative_copy_chrc
    #    number of chloroplast for each copy of nulcear genome
    # 3. Conversion rate
    
    def __init__(self, conv_rate, relative_copy_chrc, depth, n_sample = 10000):
        self.n_sample = n_sample
        self._conv_rate = conv_rate
        self._depth = depth
        self._rel_copy_chrc = relative_copy_chrc
        
        self.rsample = self.sample(conv_rate, relative_copy_chrc, depth)
        
    def sample(self, conv_rate, relative_copy_chrc, depth):
        relative_copy_chrc = relative_copy_chrc.repeat(self.n_sample,1)
        depth = depth.repeat(self.n_sample, 1)
        
        num_nulcear =  dist.Binomial(depth, 1 / (1 + relative_copy_chrc)  ).sample()
        num_non_conv = dist.Binomial(depth -  num_nulcear,  conv_rate).sample()
        return(num_nulcear + num_non_conv)
              
    def log_prob(self, x):
        # Empirical p-value
        if self.rsample.shape[1] == x.shape[0]:
            p = np.zeros(0)
            for ef_ix in np.arange(x.shape[0]):
                p = np.append(p, int((self.rsample[:,ef_ix] == x[ef_ix]).sum()) )
            p = torch.tensor(p)
        else:
            p = (self.rsample == x).sum(0)
        p[p == 0] = 1
        return( (p / self.n_sample).log() )


def get_best_conv_rate(mc_counts, mc_total, cnv_chrc):
    mc_counts = mc_counts[mc_total > 0]
    mc_total = mc_total[mc_total > 0]
    num_pos = mc_total.shape[0]
    
    outdata = {}
    outdata['conv_rates'] = torch.linspace(0.0001, 0.08, 100)
    outdata['posterior'] = np.zeros(0)
    for ef_conv in outdata['conv_rates']:
        ef_pos_ix = np.sort(np.random.randint(0, num_pos, 700))
        ef_model = CytosineModel(conv_rate=ef_conv, relative_copy_chrc=torch.tensor(int(cnv_chrc)), depth = torch.tensor(mc_total)[ef_pos_ix] )
        outdata['posterior'] = np.append(outdata['posterior'], float(ef_model.log_prob( mc_counts[ef_pos_ix] ).sum()))
    outdata['best_conv_rate'] = outdata['conv_rates'][np.argmax(outdata['posterior'])]
    return(outdata)