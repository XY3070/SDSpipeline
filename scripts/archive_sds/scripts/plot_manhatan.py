import pandas as pd

import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import kaleido

# 加载hg38注释文件
hg38_annotation_file = r"d:\pythonProject\SDS\manhattan\hg38.symbol_gene.txt"
try:
    gene_annotations = pd.read_csv(hg38_annotation_file, sep='\t', header=None, names=['CHR_hg38', 'START', 'END', 'GENE'])
    gene_annotations['CHR_hg38'] = gene_annotations['CHR_hg38'].str.replace('chr', '')
    # 过滤掉非数字染色体，并转换为整数
    gene_annotations = gene_annotations[gene_annotations['CHR_hg38'].str.isnumeric()].copy()
    gene_annotations['CHR_hg38'] = gene_annotations['CHR_hg38'].astype(int)

except FileNotFoundError:
    print(f"Error: Annotation file not found at {hg38_annotation_file}")
    gene_annotations = pd.DataFrame()
except Exception as e:
    print(f"Error loading annotation file: {e}")
    gene_annotations = pd.DataFrame()











sds_base_dir = r"D:\02 Research\2510 Phylogeny\results\raw_2511"
all_dfs = []

for i in range(1, 23): # 染色体1到22
    file_path = f"{sds_base_dir}\\chr{i}_sds_res.txt"
    try:
        # 假设SDS文件是制表符分隔的，并且有标题行
        # 请根据实际文件格式调整sep和header参数
        chr_df = pd.read_csv(file_path, sep='\t')
        # 从文件名中提取染色体号
        chr_num = int(file_path.split('chr')[1].split('_sds_res.txt')[0])
        chr_df['CHR'] = chr_num
        # 重命名列以匹配脚本预期
        chr_df.rename(columns={'ID': 'MarkerName', 'POS': 'BP'}, inplace=True)
        all_dfs.append(chr_df)
    except FileNotFoundError:
        print(f"Warning: File not found for chromosome {i}: {file_path}")
    except Exception as e:
        print(f"Error reading file for chromosome {i}: {file_path} - {e}")

if all_dfs:
    df = pd.concat(all_dfs, ignore_index=True)
else:
    raise ValueError("No SDS data files were loaded. Please check the directory and file names.")


# 2. 数据预处理
# 提取染色体号和位置（假设MarkerName格式为"chr1:1234"或"rs123"）
df['HetDf'] = 1 # 默认值，因为SDS文件中没有此列
df['Effect'] = 0 # 默认值，因为SDS文件中没有此列

# 转换数据类型（处理缺失值）
df['CHR'] = pd.to_numeric(df['CHR'], errors='coerce').astype('Int64')
df['BP'] = pd.to_numeric(df['BP'], errors='coerce').astype('Int64')
# 删除无法解析位置的标记（如非chr格式的rsID）
df = df.dropna(subset=['CHR', 'BP']).copy()

# R脚本中的DAF分箱和SDS_Z计算
# 确保DAF列存在且为数值类型
if 'DAF' not in df.columns:
    raise ValueError("DAF column not found in SDS data. Please ensure your SDS files contain a 'DAF' column.")
df['DAF'] = pd.to_numeric(df['DAF'], errors='coerce')
df = df.dropna(subset=['DAF']) # 移除DAF为NaN的行

# DAF分箱
df['daf_bin'] = pd.cut(df['DAF'], bins=np.arange(0, 1.01, 0.01), include_lowest=True, labels=False)

# 计算SDS_Z
# 确保rSDS列存在且为数值类型
if 'rSDS' not in df.columns:
    # 如果P-value列实际上是rSDS，则重命名
    if 'P-value' in df.columns:
        df.rename(columns={'P-value': 'rSDS'}, inplace=True)
    else:
        raise ValueError("rSDS or P-value column not found in SDS data. Please ensure your SDS files contain 'rSDS' or 'P-value' column.")

df['rSDS'] = pd.to_numeric(df['rSDS'], errors='coerce')
df = df.dropna(subset=['rSDS']) # 移除rSDS为NaN的行

df['SDS_Z'] = df.groupby('daf_bin')['rSDS'].transform(lambda x: (x - x.mean()) / x.std())
df['SDS_Z'] = df['SDS_Z'].fillna(0) # 填充因std为0导致的NaN

# 从SDS_Z计算P_VALUE
# 双尾P值：2 * (1 - CDF(abs(SDS_Z)))
df['P_VALUE'] = 2 * stats.norm.sf(np.abs(df['SDS_Z']))
df['P_VALUE'] = df['P_VALUE'].fillna(1) # 填充因SDS_Z为NaN导致的P_VALUE为NaN

# 使用P_VALUE计算-log10(p)
df['P_VALUE'] = df['P_VALUE'].replace(0, np.finfo(float).eps) # 将P-value为0的值替换为机器最小浮点数
df = df[df['P_VALUE'] > 0] # 过滤掉P-value小于等于0的值
df['logp'] = -np.log10(df['P_VALUE'])

# 按染色体排序
df = df.sort_values(['CHR', 'BP'])

# 3. 计算曼哈顿图x轴位置
chrom_lengths = df.groupby('CHR')['BP'].max()
df['xpos'] = df['BP'] + df['CHR'].map({c: sum(chrom_lengths.loc[:c-1]) for c in chrom_lengths.index})






def generate_manhattan_plot(df_to_plot, title, output_filename, output_formats=["html", "png"], display_threshold=None):
    df_plot = df_to_plot.copy() # Create a copy to avoid modifying the original DataFrame

    if display_threshold is not None:
        df_plot = df_plot[df_plot['P_VALUE'] <= display_threshold]
        title = f"{title} (P <= {display_threshold:.0e})" # Update title to reflect filtering
        if df_plot.empty:
            print(f"Warning: No data points to display for {title} after applying display_threshold {display_threshold}.")
            return

    fig = go.Figure()

    # 颜色方案（23条染色体）
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'] * 3

    # 添加每条染色体的数据
    for chrom in sorted(df_plot['CHR'].unique()):
        if pd.isna(chrom):
            continue
        mask = df_plot['CHR'] == chrom
        fig.add_trace(go.Scatter(
            x=df_plot[mask]['xpos'],
            y=df_plot[mask]['logp'],
            mode='markers',
            marker=dict(
                color=colors[int(chrom)-1],
                size=6,
                opacity=0.6
            ),
            name=f'Chr {int(chrom)}',
            text=df_plot[mask]['MarkerName'],
            hovertemplate=
            '<b>%{text}</b><br>' +
            'Chr: %{customdata[0]}<br>' +
            'Position: %{customdata[1]:,}<br>' +
            '-log10(p): %{y:.2f}<br>' +
            'Effect: %{customdata[2]:.3f}<extra></extra>',
            customdata=np.stack((
                df_plot[mask]['CHR'],
                df_plot[mask]['BP'],
                df_plot[mask]['Effect']
            ), axis=-1)
        ))

    # 添加显著线
    fig.add_hline(
        y=-np.log10(5e-8),
        line_dash="dash",
        line_color="red",
        annotation_text="Genome-wide significance (p=5e-8)",
        annotation_position="top right"
    )

def get_overlapping_genes(chrom, pos, gene_annotations_df):
    if gene_annotations_df.empty:
        return "No Annotation Data"
    
    # Filter annotations for the given chromosome
    chr_annotations = gene_annotations_df[gene_annotations_df['CHR_hg38'] == chrom]
    
    # Find genes where SNP position is within gene start and end
    overlapping_genes = chr_annotations[
        (chr_annotations['START'] <= pos) &
        (chr_annotations['END'] >= pos)
    ]['GENE'].tolist()
    
    if overlapping_genes:
        return ", ".join(overlapping_genes)
    else:
        return "No Gene"

    # 添加建议线（可选）
    fig.add_hline(
        y=-np.log10(1e-5),
        line_dash="dot",
        line_color="orange",
        annotation_text="Suggestive threshold (p=1e-5)",
        annotation_position="top right"
    )

    # 更新布局
    fig.update_layout(
        title=title,
        xaxis_title='Chromosome',
        yaxis_title='-log10(p-value)',
        hovermode='closest',
        showlegend=True,
        xaxis=dict(
            tickmode='array',
            tickvals=df_to_plot.groupby('CHR')['xpos'].median().values,
            ticktext=[f'{int(c)}' for c in sorted(df_to_plot['CHR'].unique()) if pd.notna(c)]
        ),
        height=600,
        width=1100
    )

    # 显示图形
    for fmt in output_formats:
        if fmt == "html":
            fig.write_html(f"{output_filename}.html")
            print(f"Manhattan plot saved as {output_filename}.html")
        elif fmt == "png":
            fig.write_image(f"{output_filename}.png", scale=2) # scale参数可以调整图片分辨率
            print(f"Manhattan plot saved as {output_filename}.png")
        else:
            print(f"Invalid output format specified: {fmt}. Skipping.")




# 定义显著性阈值
sig_threshold = 5e-8 # 显著性阈值，与R脚本保持一致
suggested_threshold = 1e-5 # 建议显著性阈值，用于图表显示

# 生成所有AF的曼哈顿图
generate_manhattan_plot(df, 'Meta-analysis Manhattan Plot (All AF)', 'manhattan_plot_all_af', output_formats=['html'], display_threshold=suggested_threshold)

# 生成AF > 1%的曼哈顿图
df_af_gt_01 = df[df['DAF'] > 0.01].copy()
generate_manhattan_plot(df_af_gt_01, 'Meta-analysis Manhattan Plot (AF > 1%)', 'manhattan_plot_af_gt_01', output_formats=['html'], display_threshold=suggested_threshold)

# 生成AF > 5%的曼哈顿图
df_af_gt_05 = df[df['DAF'] > 0.05].copy()
generate_manhattan_plot(df_af_gt_05, 'Meta-analysis Manhattan Plot (AF > 5%)', 'manhattan_plot_af_gt_05', output_formats=['html'], display_threshold=suggested_threshold)


# 5. 显著区域分析 (Significant Region Analysis)

sig_snps = df[df['P_VALUE'] < sig_threshold].copy()

if not sig_snps.empty:
    print(f"Found {len(sig_snps)} significant SNPs at threshold {sig_threshold}")
    # 按照染色体和位置排序
    sig_snps = sig_snps.sort_values(by=['CHR', 'BP']).reset_index(drop=True)

    # 添加基因注释
    sig_snps['Overlapping_Genes'] = sig_snps.apply(lambda row: get_overlapping_genes(row['CHR'], row['BP'], gene_annotations), axis=1)

    # 手动聚类逻辑
    sig_snps['cluster'] = 0
    cluster_id = 1

    if len(sig_snps) > 0:
        sig_snps.loc[0, 'cluster'] = cluster_id
        prev_chr = sig_snps.loc[0, 'CHR']
        prev_pos = sig_snps.loc[0, 'BP']

        for i in range(1, len(sig_snps)):
            curr_chr = sig_snps.loc[i, 'CHR']
            curr_pos = sig_snps.loc[i, 'BP']

            # 距离 < 200kb (200000 bp)
            if curr_chr == prev_chr and (curr_pos - prev_pos) <= 200000:
                sig_snps.loc[i, 'cluster'] = cluster_id
            else:
                cluster_id += 1
                sig_snps.loc[i, 'cluster'] = cluster_id
            prev_chr = curr_chr
            prev_pos = curr_pos

    # 汇总聚类
    final_table = sig_snps.groupby('cluster').agg(
        Region=('BP', lambda x: f"{sig_snps.loc[x.index[0], 'CHR']}:{x.min()}-{x.max()}"),
        Sig_SNVs=('MarkerName', 'count'),
        Top_SNP_ID=('P_VALUE', lambda x: sig_snps.loc[x.idxmin(), 'MarkerName']),
        Top_SNP_AF=('P_VALUE', lambda x: sig_snps.loc[x.idxmin(), 'DAF']),
        Top_SNP_SDS_Z=('P_VALUE', lambda x: round(sig_snps.loc[x.idxmin(), 'SDS_Z'], 4)),
        Top_SNP_P=('P_VALUE', 'min'),
        Genes=('Overlapping_Genes', lambda x: ", ".join(sorted(list(set(", ".join(x).split(', ')) - {"No Gene", "No Annotation Data"}))) if x.any() else "No Gene")
    ).reset_index()

    # 格式化Top_SNP_P为科学计数法
    final_table['Top_SNP_P'] = final_table['Top_SNP_P'].apply(lambda x: f'{x:.3e}')


    # 删除cluster列
    final_table = final_table.drop(columns=['cluster'])

    output_prefix = "SDS_Analysis_Result"
    summary_file = f"{output_prefix}_Significant_Regions.csv"
    final_table.to_csv(summary_file, index=False)
    print(f"Significant regions saved to {summary_file}")
    print(final_table.head())

else:
    print(f"No significant regions found at threshold {sig_threshold}")