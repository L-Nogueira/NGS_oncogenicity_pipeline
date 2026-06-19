import pandas as pd
import sys
import os
import re

def extrair_hgvs(row, coluna_fonte):
    """Função auxiliar para extrair cDNA (c.) e Proteína (p.) da string do ANNOVAR."""
    val = str(row.get(coluna_fonte, ''))
    if not val or val == '.' or val == 'nan':
        return '.', '.'
    
    # O ANNOVAR costuma separar transcritos por vírgula, pegamos o primeiro/principal
    primeiro_transcrito = val.split(',')[0]
    
    c_dot = '.'
    p_dot = '.'
    
    # Expressões regulares para capturar padrões HGVS c. e p.
    match_c = re.search(r'(c\.[^:]+)', primeiro_transcrito)
    match_p = re.search(r'(p\.[^:]+)', primeiro_transcrito)
    
    if match_c:
        c_dot = match_c.group(1)
    if match_p:
        p_dot = match_p.group(1)
        
    return c_dot, p_dot

def converter_cancervar(input_txt, output_dir, id_amostra):
    print(f"📊 Lendo e corrigindo arquivo: {os.path.basename(input_txt)}")
    
    try:
        df = pd.read_csv(input_txt, sep='\t', low_memory=False)
        
        if df.iloc[0, 1] == 'Start':
            df = df.drop(df.index[0]).reset_index(drop=True)
            
        df.columns = [c.strip() for c in df.columns]
        
        if df.columns[0].startswith('CancerVar:'):
            df.rename(columns={df.columns[0]: 'Interpretation_Details'}, inplace=True)

    except Exception as e:
        print(f"❌ Erro ao ler e processar o arquivo: {e}")
        return

    # ======================================================================
    # EXTRAÇÃO AUTOMÁTICA DE HGVS
    # ======================================================================
    print("🧬 Extraindo nomenclaturas HGVS isoladas de cDNA e Proteína...")
    
    # Tenta extrair primeiro via RefSeq (prioridade clínica), se não houver vai de Ensembl
    if 'AAChange.refGene' in df.columns:
        hgvs_dados = df.apply(lambda r: extrair_hgvs(r, 'AAChange.refGene'), axis=1)
    elif 'AAChange.ensGene' in df.columns:
        hgvs_dados = df.apply(lambda r: extrair_hgvs(r, 'AAChange.ensGene'), axis=1)
    else:
        hgvs_dados = [('.', '.')] * len(df)

    df['HGVS_cDNA'] = [h[0] for h in hgvs_dados]
    df['HGVS_Proteina'] = [h[1] for h in hgvs_dados]
    # ======================================================================

    colunas = list(df.columns)
    
    # Ordem prioritária de colunas atualizada incluindo os novos campos HGVS criados
    colunas_prioritarias = [
        'Interpretation_Details',
        'Chr', 'Start', 'End', 'Ref', 'Alt',
        'Gene.refGene', 'HGVS_cDNA', 'HGVS_Proteina',  # <-- Posicionados em destaque
        'Func.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
        'Gene.ensGene', 'AAChange.ensGene',
        'clinvar: Clinvar', 'cosmic91'
    ]
    
    colunas_para_o_inicio = [c for c in colunas_prioritarias if c in colunas]
    outras_colunas = [c for c in colunas if c not in colunas_para_o_inicio]
    df_reorganizado = df[colunas_para_o_inicio + outras_colunas]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"📁 Nova pasta criada para os relatórios: {output_dir}")

    csv_out = os.path.join(output_dir, f"{id_amostra}_relatorio_final.csv")
    xlsx_out = os.path.join(output_dir, f"{id_amostra}_relatorio_final.xlsx")

    # Salva em CSV
    df_reorganizado.to_csv(csv_out, index=False, sep=';')
    print(f"✅ CSV estruturado gerado: {csv_out}")

    # Salva em XLSX (Excel)
    with pd.ExcelWriter(xlsx_out, engine='openpyxl') as writer:
        df_reorganizado.to_excel(writer, index=False, sheet_name='Variantes_Somáticas')
    print(f"✅ Planilha Excel (.xlsx) gerada: {xlsx_out}")

if __name__ == "__main__":
    input_padrao = "/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/05_anotacao/1181_S12_L001_cancervar.output.hg38_multianno.txt.cancervar"
    dir_out_padrao = "/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/06_relatorios_finais/1181_S12_L001"
    id_padrao = "1181_S12_L001"
    
    path_in = sys.argv[1] if len(sys.argv) > 1 else input_padrao
    path_out_dir = sys.argv[2] if len(sys.argv) > 2 else dir_out_padrao
    sample_id = sys.argv[3] if len(sys.argv) > 3 else id_padrao

    if os.path.exists(path_in):
        converter_cancervar(path_in, path_out_dir, sample_id)
    else:
        print(f"❌ Arquivo de entrada não encontrado: {path_in}")
