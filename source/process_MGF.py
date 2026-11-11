import os
import re

import pandas as pd


def writeMGFfile(ouputMGFpath, originMGFpath, specname_list, mgf_keydict):
    count = 0
    flag = False
    key_spec = ""
    mgf_str = ""
    spec_name_list = []
    with open(ouputMGFpath, "a") as outMGF:
        with open(originMGFpath, "r") as oriMGF:
            for line in oriMGF:
                if line.startswith("TITLE"):
                    nline = line.strip()
                    spec_name = nline.split("=")[1]
                    key_spec = spec_name.split(".")[0]+"."+spec_name.split(".")[1]+"."+spec_name.split(".")[2]+"."+spec_name.split(".")[3]
                    if key_spec in specname_list:
                        flag = True
                        specname_list.remove(key_spec)
                        # outMGF.write("BEGIN IONS\n" + line)
                        # count += 1
                elif line.startswith("BEGIN IONS"):
                    if flag:
                        for x in mgf_keydict[key_spec]:
                            spec_name_write = key_spec + "." + x + ".dta"
                            write_str = (
                                "BEGIN IONS\nTITLE=" + spec_name_write + "\n" + mgf_str
                            )
                            spec_name_list.append(spec_name_write)
                            count += 1
                            outMGF.write(write_str)
                        mgf_str = ""
                    flag = False
                elif flag:
                    mgf_str += line

    return count, spec_name_list

def find_mgf_file(folder, key):
    pattern = re.compile(rf"^{re.escape(key)}(_[^.]+)?\.mgf$", re.IGNORECASE)
    for fname in os.listdir(folder):
        if pattern.match(fname):
            return os.path.join(folder, fname)
    return None
    
def get_dataMGF(foleder, MGFfilepath, pSite2inputfp):
    custom_dict = pd.read_csv(pSite2inputfp, sep="\t")

    res_set = set(custom_dict["spec_name"])
    res_dict = {}
    mgf_keydict = {}
    for specname in res_set:
        key_spec = specname.split(".")[0]+"."+specname.split(".")[1]+"."+specname.split(".")[2]+"."+specname.split(".")[3]
        spec_id = specname.split(".")[4]
        if key_spec in mgf_keydict.keys():
            mgf_keydict[key_spec].append(spec_id)
        else:
            mgf_keydict[key_spec] = [spec_id]

        rawid = specname.split(".")[0]
        if rawid in res_dict.keys():
            if key_spec in res_dict[rawid]:
                continue
            res_dict[rawid].append(key_spec)
        else:
            res_dict[rawid] = [key_spec]

    with open(MGFfilepath, "a+", encoding="utf-8") as test:
        test.truncate(0)
    mgf_num = 0
    specname_inmgf_list = []
    for key, value in res_dict.items():
        # file = foleder + "\\" + key + "_HCDFT.mgf"
        file = find_mgf_file(foleder, key)
        num, specname_inmgf = writeMGFfile(MGFfilepath, file, value, mgf_keydict)
        mgf_num += num
        specname_inmgf_list += specname_inmgf
        print(mgf_num)


def process_MGF(param):
    if not param.isMGFfolder:
        return
    if param.ispFind:
        get_dataMGF(param.MGFPath, param.AutoMGFPath, param.AutoCustomInputPath)
    else:
        get_dataMGF(param.MGFPath, param.AutoMGFPath, param.InputPath)
