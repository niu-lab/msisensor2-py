# coding: utf-8

import sys
import os
import xgboost


# Model prediction threshold
site_model_threshold = 0.3
msiscore_threshold = 0.25

rep_start = 1
rep_end = 26


class MSIPredict(object):
    
    def __init__(self):
        pass
    
    def normalization(self, old_list):
        sum_v = sum(old_list)
        if sum_v == 0:
            return
        for i in range(0, len(old_list)):
            old_list[i] = float(old_list[i]) / sum_v

    def load(self, file_name, site_file):
        site_dict = {}
        in_file = open(site_file, "r")
        for l in in_file:
            site_dict[l.strip()] = 1
        feature_dict = {}
        in_file_feature = open(file_name, "r")
        l = in_file_feature.readline()
        while l:
            data = l.strip().split(" ")
            loc = data[0] + "_" + data[1]
            if site_dict.get(loc, 0) == 1:
                l = in_file_feature.readline()
                l = in_file_feature.readline()
                if l[0] != "T":
                    print "distribution file error..."
                    exit()
                feature = l.strip().split(" ")[rep_start: rep_end]
                feature = map(float, feature)
                if sum(feature) < 20:
                    l = in_file_feature.readline()
                    continue
                self.normalization(feature)
                l = in_file_feature.readline()
                feature_dict[loc] = feature
            else:
                l = in_file_feature.readline()
                l = in_file_feature.readline()
                l = in_file_feature.readline()
        return feature_dict
                
    def predict(self, feature_dict, model_file):
        #predict test
        sta_num = 0
        uns_num = 0
        for site, x_test in feature_dict.items():
            bst = xgboost.Booster()
            if not os.path.exists(model_file + '/trainsites_xgb_%s.model' % site):
                #print "model missing: %s..." % site
                continue
            bst.load_model(model_file + '/trainsites_xgb_%s.model' % site)
            dtest = xgboost.DMatrix(x_test)
            y_pred = bst.predict(dtest)
            if y_pred > site_model_threshold:
                uns_num += 1
            else:
                sta_num += 1
        msiscore = float(uns_num)/(sta_num + uns_num)
        msi_status = "MSI" if float(uns_num)/(sta_num + uns_num) >= msiscore_threshold else "MSS"
        print msiscore, uns_num, sta_num
        return msi_status
    
    def sample_predict(self, file_name, site_file, model_file):
        tmp_dict = self.load(file_name, site_file)
        return self.predict(tmp_dict, model_file)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print "python msisensor2.py output_dis sites_list.txt models"
        exit()
    else:
        mst = MSIPredict()
        mst.sample_predict(sys.argv[1], sys.argv[2], sys.argv[3])





