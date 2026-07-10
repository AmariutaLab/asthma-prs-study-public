import logging

import pandas as pd
import pyfocus as pf

from sqlalchemy import create_engine


__all__ = ["import_fusion", "import_predixcan"]

# TODO: implement exporting to predixcan/fusion

COUNT = 250

def import_fusion(path, name, tissue, assay, use_ens_id, from_gencode, rsid_table, session):
    """
    Import weights from a FUSION Rdata into the FOCUS db.
    """
    log = logging.getLogger(pf.LOG)

    import re
    import os
    import warnings
    from collections import defaultdict
    import numpy as np

    try:
        import mygene
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import rpy2.robjects as robj
    except ImportError:
        log.error("Import submodule requires mygene and rpy2 to be installed.")
        raise

    log.info(f"Starting import from FUSION database {path}")
    db_ref_panel = pf.RefPanel(ref_name=name, tissue=tissue, assay=assay)
    ses = None

    load_func = robj.r['load']
    local_dir = os.path.dirname(os.path.abspath(path))
    mg = mygene.MyGeneInfo()

    fusion_db = pd.read_csv(path, delim_whitespace=True)
    if use_ens_id and from_gencode:
        genes = fusion_db.ID.apply(lambda x: re.sub("(.*)\\.\\d+", "\\1", x))
    else:
        genes = fusion_db.ID.values

    log.info("Querying mygene servers for gene annotations")
    if use_ens_id:
        results = mg.querymany(
            genes, scopes="ensembl.gene", verbose=False,
            fields="ensembl.gene,genomic_pos,symbol,ensembl.type_of_gene,alias",
            species="human"
        )
    else:
        results = mg.querymany(
            genes, scopes="symbol,alias", verbose=False,
            fields="ensembl.gene,genomic_pos,symbol,ensembl.type_of_gene,alias",
            species="human"
        )

    res_map = defaultdict(list)
    for r in results:
        res_map[r["query"]].append(r)

    dict_rsid_table = None
    if rsid_table is not None:
        log.info(f"Loading rsid table {rsid_table}")
        df = pd.read_csv(rsid_table, sep="\t")
        dict_rsid_table = {
            (str(c), p): s for c, p, s in
            zip(df["CHR"], df["POS"], df["SNP"])
        }

    DIR_IN_HEADER = "DIR" in fusion_db.columns
    count = 0
    log.info("Starting individual model conversion")

    for _, row in fusion_db.iterrows():
        wgt_dir = row.DIR if DIR_IN_HEADER else local_dir
        wgt_name, g_name, chrom, txstart, txstop = row.WGT, row.ID, row.CHR, row.P0, row.P1
        wgt_path = f"{wgt_dir}/{wgt_name}"

        # 🔥 Clear R environment to avoid object leakage
        robj.r("rm(list = ls())")
        load_func(wgt_path)

        # -------------------------
        # Gene annotation (unchanged)
        # -------------------------
        gene_info = {}
        id_dict = {}

        lookup_g_name = (
            re.sub("(.*)\\.\\d+", "\\1", g_name)
            if use_ens_id and from_gencode else g_name
        )

        for hit in res_map[lookup_g_name]:
            if "notfound" in hit or "ensembl" not in hit or "genomic_pos" not in hit:
                continue
            if not use_ens_id and hit["symbol"] != g_name and g_name not in hit.get("alias", []):
                continue

            ens = hit["ensembl"]
            pos = hit["genomic_pos"]
            if isinstance(ens, dict): ens = [ens]
            if isinstance(pos, dict): pos = [pos]

            for e in ens:
                id_dict[e["gene"]] = e["type_of_gene"]

            for p in pos:
                if not re.match("[0-9]{1,2}|X|Y", p["chr"], re.I):
                    continue
                g_id = p["ensemblgene"]
                g_type = id_dict.get(g_id)
                if not gene_info or "protein" in (g_type or ""):
                    gene_info["geneid"] = g_id
                    gene_info["type"] = g_type

        if not gene_info:
            log.warning(f"Unable to match {g_name} to Ensembl ID. Using symbol for ID")
            gene_info["geneid"] = g_name
            gene_info["type"] = None

        gene_info.update({
            "txid": None,
            "name": g_name,
            "chrom": chrom,
            "txstart": txstart,
            "txstop": txstop
        })

        # -------------------------
        # CV vs META decision
        # -------------------------
        has_cv = bool(robj.r("exists('cv.performance')"))
        wgts_raw = np.array(robj.r['wgt.matrix'])
        is_single_col = (
            wgts_raw.ndim == 1 or
            (wgts_raw.ndim == 2 and wgts_raw.shape[1] == 1)
        )

        # =========================
        # ORIGINAL IF BLOCK (UNCHANGED)
        # =========================
        if has_cv and not is_single_col:
            methods = np.array(robj.r['cv.performance'].colnames)
            types = list(robj.r['cv.performance'].rownames)

            if "rsq" not in types or "pval" not in types:
                raise ValueError(f"cv.performance missing rsq/pval for {path}")

            wgts = np.array(robj.r['wgt.matrix'])
            keep = ~np.isnan(np.std(wgts, axis=0))
            wgts = wgts.T[keep].T
            methods = methods[keep]

            values = np.array(robj.r['cv.performance'])
            if values.shape[0] > values.shape[1]:
                values = values.T
            values = values.T[keep].T

            rsq_idx = types.index("rsq")
            pval_idx = types.index("pval")

            r2 = -100
            pval = 1
            r2idx = 0
            method = None
            top1_idx = -1

            for i, v in enumerate(values[rsq_idx]):
                if methods[i] == "top1":
                    top1_idx = i
                    continue
                if v > r2:
                    r2 = v
                    method = methods[i]
                    r2idx = i
                    pval = values[pval_idx, i]

            if method is None:
                method = "top1"
                r2idx = top1_idx
                r2 = values[rsq_idx][top1_idx]
                pval = values[pval_idx][top1_idx]

            wgts = wgts.T[r2idx]
            attrs = {"cv.R2": r2, "cv.R2.pval": pval}

        # =========================
        # ELSE: META MODEL ONLY
        # =========================
        else:
            wgts = wgts_raw[:, 0] if wgts_raw.ndim == 2 else wgts_raw
            keep_nan = ~np.isnan(wgts)
            wgts = wgts[keep_nan]
            method = "meta"
            attrs = {"cv.R2": None, "cv.R2.pval": None}

        # -------------------------
        # SNP handling (unchanged)
        # -------------------------
        snps = robj.r['snps']
        snp_info = pd.DataFrame({
            "snp": list(snps[1]),
            "chrom": [str(chrom) for chrom in snps[0]],
            "pos": list(snps[3]),
            "a1": list(snps[4]),
            "a0": list(snps[5])
        })

        if dict_rsid_table is not None:
            snp_info["snp"] = snp_info.apply(
                lambda r: dict_rsid_table.get((r.chrom, r.pos)), axis=1
            )
            keep = snp_info.snp.notnull().values
            wgts = wgts[keep]
            snp_info = snp_info[keep]

        keep = ~np.isclose(wgts, 0)
        wgts = wgts[keep]
        snp_info = snp_info[keep]

        model = pf.build_model(
            gene_info, snp_info, db_ref_panel, wgts, ses, attrs, method
        )

        session.add(model)
        try:
            session.commit()
        except Exception:
            session.rollback()
            raise

        count += 1
        if count % COUNT == 0:
            log.info(f"Committed {COUNT} models to db")

    if count % COUNT:
        log.info(f"Committed {count % COUNT} models to db")

    log.info(f"Finished import from FUSION database {path}")
    return


def export_fusion(path, session):
    log = logging.getLogger(pf.LOG)
    raise NotImplementedError("export_fusion not implemented!")
    return


def import_predixcan(path, name, tissue, assay, method, session):
    """
    Import weights from a PrediXcan db into the FOCUS db.

    :param path:  string path to the PrediXcan db
    :param name: str name of the reference panel
    :param tissue: str name of the tissue
    :param assay: technology assay to measure abundance
    :param method: the prediction model used to fit the data
    :param session: sqlalchemy.Session object for the FOCUS db

    :return:  None
    """
    log = logging.getLogger(pf.LOG)

    import re
    import os
    import numpy as np

    from collections import defaultdict
    try:
        import mygene
    except ImportError:
        log.error("Import submodule requires mygene and rpy2 to be installed.")
        raise

    if not os.path.isfile(path):
        raise ValueError(f"Cannot find database {path}")
    log.info(f"Starting import from PrediXcan database {path}")
    pred_engine = create_engine(f"sqlite:///{path}")

    weights = pd.read_sql_table('weights', pred_engine)
    extra = pd.read_sql_table('extra', pred_engine)

    def gencode2ensmble(x):
        idx = x.rfind(".")
        return x if idx == -1 else x[:idx]

    # get unique genes
    genes = weights.gene.unique()
    genes = [gencode2ensmble(g) for g in genes]

    log.info("Querying mygene servers for gene annotations")
    mg = mygene.MyGeneInfo()
    results = mg.querymany(genes, scopes='ensembl.gene', verbose=False,
                           fields="genomic_pos_hg19,symbol,alias", species="human")

    res_map = defaultdict(list)
    for result in results:
        res_map[result["query"]].append(result)

    db_ref_panel = pf.RefPanel(ref_name=name, tissue=tissue, assay=assay)
    ses = None

    count = 0
    log.info("Starting individual model conversion")
    for gid, gene in weights.groupby("gene"):
        log.debug(f"Importing gene model {gid}")
        gene_extra = extra.loc[extra.gene == gid]

        chrom = gene.varID.values[0].split("_")[0]  # grab chromosome from varID
        pos = gene.varID.map(lambda x: int(x.split("_")[1])).values  # grab basepair pos
        txstart = txstop = np.median(pos)

        g_id = gene_extra.gene.values[0]
        g_name = gene_extra.genename.values[0]
        query_id = gencode2ensmble(g_id)

        for hit in res_map[query_id]:
            if "notfound" in hit:
                continue

            if hit["symbol"] != g_name and "alias" in hit and g_name not in hit["alias"]:
                continue

            if "genomic_pos_hg19" not in hit:
                continue

            gpos = hit["genomic_pos_hg19"]
            if type(gpos) is dict:
                gpos = [gpos]

            for entry in gpos:
                # skip non-primary assembles. they have weird CHR entries e.g., CHR_HSCHR1_1_CTG3
                if not re.match("[0-9]{1,2}|X|Y", entry["chr"], re.IGNORECASE):
                    continue

                txstart = entry['start']
                txstop = entry['end']
                break

            if txstart is not None:
                # we want to use standardized Ensembl identifiers; not GENCODE modified ones...
                g_id = query_id
                break

        gene_info = dict()
        gene_info["geneid"] = g_id
        gene_info["txid"] = None
        gene_info["name"] = g_name
        gene_info["type"] = gene_extra.gene_type.values[0]
        gene_info["chrom"] = chrom
        gene_info["txstart"] = txstart
        gene_info["txstop"] = txstop

        snp_info = pd.DataFrame({"snp": gene.rsid.values,
                                "chrom": [chrom] * len(gene),
                                "pos": pos,
                                "a1": gene.eff_allele.values,
                                "a0": gene.ref_allele.values})

        wgts = gene.weight.values

        attrs = dict()
        # predixcan enet uses cv_R2_avg and nested_cv_fisher_pval
        # mashr uses pred.perf.R2 and pred.perf.pval
        if "cv_R2_avg" in gene_extra:
            attrs["cv.R2"] = gene_extra["cv_R2_avg"].values[0]
            attrs["cv.R2.pval"] = gene_extra["nested_cv_fisher_pval"].values[0]
        elif "pred.perf.R2" in gene_extra:
            attrs["cv.R2"] = gene_extra["pred.perf.R2"].values[0]
            attrs["cv.R2.pval"] = gene_extra["pred.perf.pval"].values[0]

        # build model
        model = pf.build_model(gene_info, snp_info, db_ref_panel, wgts, ses, attrs, method)
        session.add(model)
        try:
            session.commit()
        except Exception as comm_err:
            session.rollback()
            raise Exception("Failed committing to db")

        count += 1
        if count % COUNT == 0:
            log.info(f"Committed {COUNT} models to db")

    if count % COUNT != 0:
        log.info(f"Committed {count % COUNT} models to db")


    log.info(f"Finished import from PrediXcan database {path}")
    return


def export_predixcan(path, session):
    log = logging.getLogger(pf.LOG)
    raise NotImplementedError("export_predixcan not implemented!")
    return
