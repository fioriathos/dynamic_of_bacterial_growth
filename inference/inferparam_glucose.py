import minimize_lengths as mn
import recursive_lengths as rl
import pandas as pd
import numpy as np
import sys
import pickle
from scipy.stats import linregress 
#file_to_analy ='/Users/fiori/DoubleAdderArticle/PreProcessed/20180706_GW296_glycerol37_r095_div.csv'
#file_to_analy\
#='/Users/fiori/DoubleAdderArticle/PreProcessed/20180709_GW296_glucose8aa_r095_div.csv'
file_to_analy='/Users/fiori/DoubleAdderArticle/PreProcessed/20180711_GW296_glucose_r095_div.csv'
#sys.argv[1]
step = int(sys.argv[1])
#step=8
# minimal cell length
fil = 15#int(sys.argv[3])
#fil=24
leng = 'length_box_um'#sys.argv[4]
dglu = pd.read_csv(file_to_analy)
dglu = rl.give_good_structure(dglu)
dft = rl.give_unique_dataset(dglu,step,fil)
_,in_dic = rl.build_data_strucutre(dft,leng,1)
boundary = [(1e-10,None),(1e-5,None),(1e-10,None),(1e-10,None),(1e-10,None)]
m,g,s,e,a =[in_dic['s'][1,0],0.01,2.3630e-07,in_dic['sm2'],in_dic['sd2']]
r = np.random.rand(5)*np.random.choice((-1,1),5)*0.5 #50% max var 
mod=mn.minimize_lengths(free={'sl2':s+s*r[0],'gamma':g+g*r[1],'sm2':e+e*r[2],'mlam':m+m*r[3],'sd2':a},fixed={},boundary=boundary)
bestpar =\
    mod.minimize(in_dic=in_dic,numerical=False,fun=rl.grad_obj_wrap)
mb = mod.errorbars(in_dic)
mb['log_lik']= bestpar['log_lik']
mb['step']=step
mb['filt']=fil
f = open("{}.pkl".format('gluc_pap{}'.format(step)),"wb")
pickle.dump(mb,f)
f.close()
#mbp = mb['param']
#mbe = mb['error']
#print('log_lik','ml','gamma','sl2','sm2','sd2','eml','egamma','esl2','esm2','esd2','step','medium','file','length','min_cell_len','discard_top','end_type','filt_r')
#print(bestpar['log_lik'],mbp['mlam'],mbp['gamma'],mbp['sl2'],mbp['sm2'],mbp['sd2'],mbe['mlam'],mbe['gamma'],mbe['sl2'],mbe['sm2'],mbe['sd2'],step,dft.medium.iloc[0],\
#          file_to_analy,leng,fil,'1','div',str(filt_r))
    #print(bestpar['log_lik'],mod.correct_scaling(in_dic),mod.errorbars(in_dic))


