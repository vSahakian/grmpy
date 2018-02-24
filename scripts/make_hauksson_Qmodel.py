############################################################################
######################  Make Q model, save to file #########################
############################################################################

# VJS 2/2018

import dread as dr
import cdefs as cdf
import cPickle as pickle

######### PATHS ###########

## Input:
hauksson_text_path_Qp = '/Users/vsahakian/anza/data/vm/Hauksson/vl_sc_q06d_hhz_Qp.plt'
hauksson_text_path_Qs = '/Users/vsahakian/anza/data/vm/Hauksson/vl_sc_q08d_hhe_Qs.plt'

## Output:
hauksson_model_path_Qp = '/Users/vsahakian/anza/data/vm/Hauksson/Qp.pckl'
hauksson_model_path_Qs = '/Users/vsahakian/anza/data/vm/Hauksson/Qs.pckl'

#Read in data:
x_Qp,y_Qp,z_Qp,nx_Qp,ny_Qp,nz_Qp,model_Qp=dr.read_hauksson_file(hauksson_text_path_Qp)
x_Qs,y_Qs,z_Qs,nx_Qs,ny_Qs,nz_Qs,model_Qs=dr.read_hauksson_file(hauksson_text_path_Qs)

#Make object:
materialobject_Qp=cdf.material_model(x_Qp,y_Qp,z_Qp,nx_Qp,ny_Qp,nz_Qp,model_Qp)
materialobject_Qs=cdf.material_model(x_Qs,y_Qs,z_Qs,nx_Qs,ny_Qs,nz_Qs,model_Qs)


#Save it to a file:
mfile=open(hauksson_model_path_Qp,'w')
pickle.dump(materialobject_Qp,mfile)
mfile.close()

#Save it to a file:
mfile=open(hauksson_model_path_Qs,'w')
pickle.dump(materialobject_Qs,mfile)
mfile.close()

