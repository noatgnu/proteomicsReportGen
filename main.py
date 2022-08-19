import io
import os
import re
import pandas as pd
from plotnine import ggplot, aes, geom_point, ggsave, geom_col, position_dodge2, theme_minimal, geom_violin, \
    geom_dotplot, coord_cartesian
from uniprotparser.betaparser import UniprotParser, UniprotSequence
import matplotlib.pyplot as plt
import seaborn as sns

os.add_dll_directory(r"C:\Program Files\GTK2-Runtime Win64\bin")

from weasyprint import HTML, CSS


differential_analysis_file = r"C:\Users\Toan Phung\Downloads\different.txt"
raw_file = r"C:\Users\Toan Phung\Downloads\raw.txt"
output_file = "test.pdf"
columns = "accession,id,gene_names,protein_name,organism_name,organism_id,length,xref_refseq,go_id,go_p,go_c,go_f,cc_subcellular_location,ft_topo_dom,ft_carbohyd,mass,cc_mass_spectrometry,sequence,ft_var_seq,cc_alternative_products,cc_function,ft_domain,xref_string"


protein_list = []
primary_id_column_differential = "Majority protein IDs"
primary_id_column_raw = "T: Majority protein IDs"

sample_column_raw = ["R1441G_HG.1", "R1441G_HG.2", "R1441G_HG.3", "R1441G_mock.1", "R1441G_mock.2", "R1441G_mock.3"]

pvalue_column = "#NAME?"
difference_column = "Difference"

domain_len = re.compile("(\d*)\.\.(\d*)")
domain_name = re.compile('note="(\.+)"')

bootstrap_link = "https://cdn.jsdelivr.net/npm/bootstrap@5.2.0/dist/css/bootstrap.min.css"

def make_html(uniprot_data, filename):


    html = HTML(string=
                f"""
    <!DOCTYPE html>
    <html>

      <head>
      
        <title>PDF Generation with Python and WeasyPrint</title>
        <link href="file:///C:/Users/Toan Phung/Documents/GitHub/proteomicsReportGen/style.css" rel="stylesheet" type="text/css">
        <link href="{bootstrap_link}" rel="stylesheet">
      </head>
      <body>
      <div class="container-fluid">
        <h1 class="text-primary">{uniprot_data["Gene Names"]}</h1>
        <div class="card">
            <div class="card-body">
                {uniprot_data["Protein names"]}
            </div>
        </div>
        <div class="card">
            <div class="card-body">
                <div class="row">
                    <div class="col-6">
                        <img style="max-width: 300px;" src="file:///C:/Users/Toan Phung/Documents/GitHub/proteomicsReportGen/{filename+'violin.png'}">
                    </div>
                </div>
                
            </div>
        </div>
        <div class="card">
            <div class="card-body">
                <img style="max-height: 300px;" src="file:///C:/Users/Toan Phung/Documents/GitHub/proteomicsReportGen/{filename+'volcano.png'}">
            </div>
        </div>
        <div class="card">
            <div class="card-body">
                <img style="max-height: 300px;" src="file:///C:/Users/Toan Phung/Documents/GitHub/proteomicsReportGen/{filename+'barchart.png'}">
            </div>
        </div>
      </div>
      </body>
    </html>
                        """)

    return html


def draw_volcano(background_df, filename, selected_df):
    fig = (ggplot() +
           geom_point(background_df, aes(x=difference_column, y=pvalue_column, fill="Category"),
                      alpha=0.1) +
           geom_point(selected_df,
                      aes(x=difference_column, y=pvalue_column, fill="Category"), alpha=1) + theme_minimal()
           )
    ggsave(fig, filename + "volcano.png")

def draw_barchart(df, filename):
    #fig = (ggplot(df, aes(x="Category", y="value", fill="Category")) + geom_col(position=position_dodge2()) + theme_minimal())
    #ggsave(fig, filename + "barchart.png")
    plt.clf()
    fig, ax= plt.subplots()
    value_count = df["Category"].value_counts()
    sns.catplot(data=df, x="Sample", y="value", hue="Category", dodge=False, kind="bar", ax=ax)
    plt.xticks(rotation=90)
    plt.tight_layout()
    tick = ax.get_xticklabels()
    labels = []
    current_pos = 0
    for n, i in enumerate(tick):

        label = i.get_text()
        temp = df[df["Sample"]==label].iloc[0]
        if value_count[temp]/2 == n+1:

        position_within = (n - current_pos) / value_count[temp]


    plt.savefig(filename + "barchart.png")

def draw_violin(df, filename):
    plt.clf()
    sns.violinplot(data=df, x="Category", y="value", hue="Category")
    sns.stripplot(data=df, x="Category", y="value", color=".3")
    #fig = (ggplot(df, aes(x="Category", y="value", fill="Category")) + geom_violin(trim=False) + coord_cartesian(expand=True))
    #ggsave(fig, filename + "violin.png")
    plt.savefig(filename + "violin.png")

if __name__ == "__main__":
    differential_df = pd.read_csv(differential_analysis_file, sep="\t")
    raw_df = pd.read_csv(raw_file, sep="\t")
    combined_df = differential_df.merge(raw_df, left_on=primary_id_column_differential, right_on=primary_id_column_raw)
    acc_set = set()
    if not os.path.exists("temp.txt"):
        for a in combined_df[primary_id_column_differential].str.split(";").explode().str.split("-"):
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
            if pd.isnull(dict_map[r[primary_id_column_differential]]["Gene Names"]):
                df_selected["Category"] = df_selected[primary_id_column_differential]
            else:
                df_selected["Category"] = dict_map[r[primary_id_column_differential]]["Gene Names"]
            df_background = combined_df[combined_df[primary_id_column_differential] != r[primary_id_column_differential]]
            df_background["Category"] = "Background"

            draw_volcano(df_background, r[primary_id_column_differential], df_selected)

            df_raw_selected = df_selected[[primary_id_column_differential] + sample_column_raw]
            df_raw = df_raw_selected.melt(id_vars=primary_id_column_differential, value_vars=sample_column_raw, var_name="Sample")
            df_raw[["Category", "Replicate"]] = df_raw["Sample"].str.split(".", expand=True)
            df_raw["Category"] = pd.Categorical(df_raw["Category"])
            draw_barchart(df_raw, r[primary_id_column_differential])
            draw_violin(df_raw, r[primary_id_column_differential])
            htmls.append(make_html(
                dict_map[r[primary_id_column_differential]],
                r[primary_id_column_differential]
                     ).render(stylesheets=[CSS(bootstrap_link)]))
    pdf_pages = []
    for h in htmls:
        for p in h.pages:
            pdf_pages.append(p)

    output_pdf = htmls[0].copy(pdf_pages).write_pdf(target=output_file)



