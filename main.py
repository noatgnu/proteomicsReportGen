import io
import os
import re
import pandas as pd
from plotnine import ggplot, aes, geom_point, ggsave
from uniprotparser.betaparser import UniprotParser, UniprotSequence
import matplotlib.pyplot as plt
from weasyprint import HTML, CSS

os.add_dll_directory(r"C:\Program Files\GTK3-Runtime Win64\bin")



differential_analysis_file = r"C:\Users\Toan Phung\Downloads\For_Curtain_SGK3-TP.txt"
raw_file = r"C:\Users\Toan Phung\Downloads\For_Curtain_Raw_SGK3_TP_Screen2.txt"
output_file = "test.pdf"
columns = "accession,id,gene_names,protein_name,organism_name,organism_id,length,xref_refseq,go_id,go_p,go_c,go_f,cc_subcellular_location,ft_topo_dom,ft_carbohyd,mass,cc_mass_spectrometry,sequence,ft_var_seq,cc_alternative_products,cc_function,ft_domain,xref_string"


protein_list = []
primary_id_column_differential = "Majority.protein.IDs"
primary_id_column_raw = "T: Majority protein IDs"

sample_column_raw = ["DMSO.01", "DMSO.02", "DMSO.03", "MK-2206.01", "MK-2206.02", "MK-2206.03", "14H.01", "14H.02", "14H.03"]

pvalue_column = "p-Value(-Log10): DMSO/14H"
difference_column = "Difference (Log2): DMSO/14H"

domain_len = re.compile("(\d*)\.\.(\d*)")
domain_name = re.compile('note="(\.+)"')

bootstrap_link = "https://cdn.jsdelivr.net/npm/bootstrap@5.2.0/dist/css/bootstrap.min.css"

def make_html(selected_df, background_df, uniprot_data, filename):
    fig = (ggplot() +
           geom_point(background_df, aes(x=difference_column, y=pvalue_column, fill="Category"),
                      alpha=0.1) +
           geom_point(selected_df,
                      aes(x=difference_column, y=pvalue_column, fill="Category"), alpha=1)
           )
    ggsave(fig, filename+".svg")
    html = HTML(string=
                f"""
    <!DOCTYPE html>
    <html>

      <head>
        <title>PDF Generation with Python and WeasyPrint</title>
        <link href="{bootstrap_link}" rel="stylesheet" />
      </head>
      <body>
      <div class="container">
        <p class="display-6">{uniprot_data["Gene Names"]}</p>
        <div class="card">
            <div class="card-body">
                <img src="file:///D:/PycharmProjects/proteomicsReportGen/{filename+'.svg'}">
            </div>
        </div>
      </div>
      </body>
    </html>
                        """)

    return html


if __name__ == "__main__":
    differential_df = pd.read_csv(differential_analysis_file, sep="\t")
    raw_df = pd.read_csv(raw_file, sep="\t")
    combined_df = differential_df.merge(raw_df, left_on=primary_id_column_differential, right_on=primary_id_column_raw)
    acc_set = set()
    if not os.path.exists("temp.txt"):
        for a in combined_df["Majority.protein.IDs"].str.split(";").explode().str.split("-"):
            acc_set.add(a[0])
        parser = UniprotParser(columns=columns)
        df = []
        for p in parser.parse(list(acc_set)):
            df.append(pd.read_csv(io.StringIO(p), sep="\t"))
        if len(df) > 1:
            df = pd.concat(df, ignore_index=True)
        else:
            df = df[0]
        df.to_csv("temp.txt", sep="\t", index=False)
    else:
        df = pd.read_csv("temp.txt", sep="\t")

    for i, r in df.iterrows():
        if pd.notnull(r["Domain [FT]"]):
            result = []
            dm = r["Domain [FT]"].split(";")
            domain = {"name": "", "start": 0, "end": 0}
            for d in dm:
                if d.startswith("DOMAIN"):
                    if domain["name"] != "":
                        result.append(domain.copy())
                    domain = {"name": "", "start": 0, "end": 0}
                    match = domain_len.search(d)
                    if match:
                        if match.group(1):
                            domain["start"] = int(match.group(1))
                        if match.group(2):
                            domain["end"] = int(match.group(2))
                elif d.startswith("/note="):
                    match = domain_name.search(d)
                    if match:
                        domain["name"] = match.group(1)
            result.append(domain)
            df.at[i, "domain_dict"] = result
    dict_map = {}

    htmls = []

    if len(protein_list) == 0:

        for i, r in combined_df.iterrows():
            uniprot_data = {}
            for p in r[primary_id_column_differential].split(";"):
                accession = UniprotSequence(p, parse_acc=True)
                temp = df[df["From"] == accession.accession]
                if len(temp) > 0:
                    uniprot_data = temp.iloc[0]
                    break
            dict_map[r[primary_id_column_differential]] = uniprot_data
        selected = combined_df.head(50)
        for i, r in selected.iterrows():
            df_selected = selected[selected[primary_id_column_differential]==r[primary_id_column_differential]]
            df_selected["Category"] = dict_map[r[primary_id_column_differential]]["Gene Names"]
            df_background = combined_df[combined_df[primary_id_column_differential]!=r[primary_id_column_differential]]
            df_background["Category"] = "Background"
            htmls.append(make_html(
                df_selected[[primary_id_column_differential, pvalue_column, difference_column, "Category"]],
                df_background[[primary_id_column_differential, pvalue_column, difference_column, "Category"]],
                dict_map[r[primary_id_column_differential]],
                r[primary_id_column_differential]
                     ).render(stylesheets=[CSS(bootstrap_link)]))
    pdf_pages = []
    for h in htmls:
        for p in h.pages:
            pdf_pages.append(p)

    output_pdf = htmls[0].copy(pdf_pages).write_pdf(target=output_file)



