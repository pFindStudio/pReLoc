class PARAM:
    def __init__(self, paramFilePath):
        self.paramFilePath = paramFilePath
        self.CparamFilePath = paramFilePath
        self.ispFindOpen = False
        self.ispChem = False
        self.pChemmodnum = 1
        self.InputSummary = ""
        self.ispPred = False
        self.pFindResfp = ""
        self.rawfolder = ""
        self.ftmodelfolder = ""
        self.ftmodelname = ""
        self.instrument = 2
        self.nce = 28
        self.rtparam = "0_100"
        self.ms2param = "0_28_28_28_1_6"
        self.score_model = "XGBOOST"
        self.cl_device = "auto"
        self.cl_batch_size = 1
        self.cl_model_path = "best_model.pth"
        self.cl_esm_path = "esm2_t33_650M_UR50D.pt"
        (
            self.isMGFfolder,
            self.MGFPath,
            self.ispFind,
            self.ispFindOpen,
            self.InputPath,
            self.ispChem,
            self.pChemmodnum,
            self.InputSummary,
            self.score_isSite,
            self.score_isMod,
            self.ispPred,
            self.pFindResfp,
            self.rawfolder,
            self.ftmodelfolder,
            self.ftmodelname,
            self.instrument,
            self.nce,
            self.rtparam,
            self.ms2param,
            self.isFLRest,
            self.FLRthresh,
            self.ModName,
            self.scoreModName,
            self.AutoCustomInputPath,
            self.DirectResPath,
            self.CandidatesPath,
            self.pSiteResPath,
            self.TimePath,
            self.CparamFilePath,
            self.AutoMGFPath,
        ) = self.read_param()

    def read_param(self):
        scoreModName = []
        ModName = []

        InputPathline = False
        modsline = False
        scoreline = False
        scoreModline = False
        Resultline = False
        pPredline = False
        FLRline = False

        with open(self.paramFilePath, "r", encoding="utf-8-sig") as paramfile:
            for line in paramfile.readlines():
                line = line.strip()
                if line.startswith("[path_info]"):
                    InputPathline = True
                    modsline = False
                    scoreline = False
                    scoreModline = False
                    Resultline = False
                    pPredline = False
                    FLRline = False
                    continue
                elif (
                    line.startswith("[mods]")
                    or line.startswith("[NTerm_mods]")
                    or line.startswith("[CTerm_mods]")
                ):
                    InputPathline = False
                    modsline = True
                    scoreline = False
                    scoreModline = False
                    Resultline = False
                    pPredline = False
                    FLRline = False
                    continue
                elif line.startswith("[score]"):
                    InputPathline = False
                    modsline = False
                    scoreline = True
                    scoreModline = False
                    Resultline = False
                    pPredline = False
                    FLRline = False
                    continue
                elif line.startswith("[scoreModName]"):
                    InputPathline = False
                    modsline = False
                    scoreline = False
                    scoreModline = True
                    Resultline = False
                    pPredline = False
                    FLRline = False
                    continue
                elif line.startswith("[output_info]"):
                    InputPathline = False
                    modsline = False
                    scoreline = False
                    scoreModline = False
                    Resultline = True
                    pPredline = False
                    FLRline = False
                    continue
                elif line.startswith("[pPred]"):
                    InputPathline = False
                    modsline = False
                    scoreline = False
                    scoreModline = False
                    Resultline = False
                    pPredline = True
                    FLRline = False
                    continue
                elif line.startswith("[FLR]"):
                    InputPathline = False
                    modsline = False
                    scoreline = False
                    scoreModline = False
                    Resultline = False
                    pPredline = False
                    FLRline = True
                    continue
                elif line.startswith("[") or line == "" or line.startswith("#"):
                    continue

                if InputPathline:
                    if line.startswith("MGFPath"):
                        self.MGFPath = line.split("=")[1]
                    elif line.startswith("InputSummary"):
                        self.InputSummary = line.split("=")[1]
                    elif line.startswith("Input"):
                        self.InputPath = line.split("=")[1]
                    elif line.startswith("isMGFfolder"):
                        self.isMGFfolder = line.split("=")[1].lower() == "true"
                    elif line.startswith("ispfindOpen"):
                        self.ispFindOpen = line.split("=")[1].lower() == "true"
                    elif line.startswith("ispfind"):
                        self.ispFind = line.split("=")[1].lower() == "true"
                    elif line.startswith("pChemmodnum"):
                        self.pChemmodnum = int(line.split("=")[1])
                    elif line.startswith("ispChem"):
                        self.ispChem = line.split("=")[1].lower() == "true"

                    else:
                        continue

                if modsline:
                    modname = line.split("=")[1]
                    ModName.append(modname)

                if scoreline:
                    key, _, value = line.partition("=")
                    key, value = key.strip().lower(), value.strip()
                    if key == "model":
                        self.score_model = value.upper()
                    elif key == "cldevice":
                        self.cl_device = value.lower()
                    elif key == "clbatchsize":
                        self.cl_batch_size = int(value)
                    elif key == "clmodelpath":
                        self.cl_model_path = value
                    elif key == "clesmpath":
                        self.cl_esm_path = value
                    if line.startswith("isSite"):
                        self.score_isSite = line.split("=")[1].lower() == "true"
                    elif line.startswith("isMod"):
                        self.score_isMod = line.split("=")[1].lower() == "true"
                    else:
                        continue

                if scoreModline:
                    modname = line.split("=")[1]
                    scoreModName.append(modname)

                if pPredline:
                    if line.startswith("pFindResfp"):
                        self.pFindResfp = line.split("=")[1]
                    elif line.startswith("rawfolder"):
                        self.rawfolder = line.split("=")[1]
                    elif line.startswith("modelfolder"):
                        self.ftmodelfolder = line.split("=")[1]
                    elif line.startswith("modelname"):
                        self.ftmodelname = line.split("=")[1]
                    elif line.startswith("instrument"):
                        self.instrument = line.split("=")[1]
                    elif line.startswith("nce"):
                        self.nce = line.split("=")[1]
                    elif line.startswith("rtparam"):
                        self.rtparam = line.split("=")[1]
                    elif line.startswith("ms2param"):
                        self.ms2param = line.split("=")[1]
                    elif line.startswith("ispPred"):
                        self.ispPred = line.split("=")[1].lower() == "true"
                    else:
                        continue

                if FLRline:
                    if line.startswith("isFLR"):
                        self.isFLRest = int(line.split("=")[1])
                    elif line.startswith("FLR"):
                        self.FLRthresh = float(line.split("=")[1])
                    else:
                        continue

                if Resultline:
                    if line.startswith("Output"):
                        OutputfolderPath = line.split("=")[1]

                        self.AutoCustomInputPath = (
                            OutputfolderPath + "\\AutoCustomInput.tsv"
                        )
                        self.DirectResPath = OutputfolderPath + "\\DirectRes.tsv"
                        self.CandidatesPath = OutputfolderPath + "\\Candidates.ld"
                        self.pSiteResPath = OutputfolderPath + "\\pReLocRes.tsv"
                        self.TimePath = OutputfolderPath + "\\Time.txt"
                        self.AutoMGFPath = OutputfolderPath + "\\dataMGF.mgf"
                        if self.isFLRest == 2 or self.ispChem:
                            self.CparamFilePath = OutputfolderPath + "\\paramtmp.txt"

        if self.score_model not in {"XGBOOST", "CL"}:
            raise ValueError("model must be XGBoost or CL")
        if self.cl_batch_size < 1:
            raise ValueError("clBatchSize must be at least 1")

        self.scoreModName = scoreModName
        self.ModName = ModName

        return (
            self.isMGFfolder,
            self.MGFPath,
            self.ispFind,
            self.ispFindOpen,
            self.InputPath,
            self.ispChem,
            self.pChemmodnum,
            self.InputSummary,
            self.score_isSite,
            self.score_isMod,
            self.ispPred,
            self.pFindResfp,
            self.rawfolder,
            self.ftmodelfolder,
            self.ftmodelname,
            self.instrument,
            self.nce,
            self.rtparam,
            self.ms2param,
            self.isFLRest,
            self.FLRthresh,
            self.ModName,
            self.scoreModName,
            self.AutoCustomInputPath,
            self.DirectResPath,
            self.CandidatesPath,
            self.pSiteResPath,
            self.TimePath,
            self.CparamFilePath,
            self.AutoMGFPath,
        )
