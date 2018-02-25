############################################################################
######################  Make Q model, save to file #########################
############################################################################

# VJS 2/2018

import dread as dr
import cdefs as cdf
import cPickle as pickle
import numpy as np

######### PATHS ###########

## Input:
hauksson_text_path_Qp = '/Users/vsahakian/anza/data/vm/Hauksson/vl_sc_q06d_hhz_Qp.plt'
hauksson_text_path_Qs = '/Users/vsahakian/anza/data/vm/Hauksson/vl_sc_q08d_hhe_Qs.plt'

## Output:
hauksson_model_path_Qp = '/Users/vsahakian/anza/data/vm/Hauksson/Qp.pckl'
hauksson_model_path_Qs = '/Users/vsahakian/anza/data/vm/Hauksson/Qs.pckl'
hauksson_model_path_Qp_simple = '/Users/vsahakian/anza/data/vm/Hauksson/Qp_simple.txt'
hauksson_model_path_Qs_simple = '/Users/vsahakian/anza/data/vm/Hauksson/Qs_simple.txt'

#Read in data:
x_Qp,y_Qp,z_Qp,nx_Qp,ny_Qp,nz_Qp,model_Qp,nodex_Qp, nodey_Qp, nodez_Qp, nodeQ_Qp = dr.read_hauksson_file(hauksson_text_path_Qp,provide_simple=True)
x_Qs,y_Qs,z_Qs,nx_Qs,ny_Qs,nz_Qs,model_Qs,nodex_Qs, nodey_Qs, nodez_Qs, nodeQ_Qs = dr.read_hauksson_file(hauksson_text_path_Qs,provide_simple=True)

#Make object:
materialobject_Qp=cdf.material_model(x_Qp,y_Qp,z_Qp,nx_Qp,ny_Qp,nz_Qp,model_Qp)
materialobject_Qs=cdf.material_model(x_Qs,y_Qs,z_Qs,nx_Qs,ny_Qs,nz_Qs,model_Qs)

# Save simple data to an array for file:
simple_model_Qp = np.c_[nodex_Qp, nodey_Qp, nodez_Qp, nodeQ_Qp]
simple_model_Qs = np.c_[nodex_Qs, nodey_Qs, nodez_Qs, nodeQ_Qs]


#Save it to a file:
mfile=open(hauksson_model_path_Qp,'w')
pickle.dump(materialobject_Qp,mfile)
mfile.close()

#Save it to a file:
mfile=open(hauksson_model_path_Qs,'w')
pickle.dump(materialobject_Qs,mfile)
mfile.close()

np.savetxt(hauksson_model_path_Qp_simple,simple_model_Qp,fmt='%.6f\t%.6f\t%.2f\t%.2f')
np.savetxt(hauksson_model_path_Qs_simple,simple_model_Qs,fmt='%.6f\t%.6f\t%.2f\t%.2f')

