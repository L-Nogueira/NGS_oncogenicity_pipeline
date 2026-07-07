import pandas as pd
import sys
import os
import re
import glob

def extrair_hgvs(row, coluna_fonte):
    """Extrai cDNA (c.) e Proteína (p.) apenas do primeiro transcrito (canônico)."""
    val = str(row.get(coluna_fonte, ''))
    if not val or val == '.' or val == 'nan':
        return '.', '.'
    
    primeiro_transcrito = val.split(',')[0]
    c_dot = '.'
    p_dot = '.'
    
    match_c = re.search(r'(c\.[^:]+)', primeiro_transcrito)
    match_p = re.search(r'(p\.[^:]+)', primeiro_transcrito)
    
    if match_c:
        c_dot = match_c.group(1)
    if match_p:
        p_dot = match_p.group(1)
        
    return c_dot, p_dot

def extrair_vaf_e_depth_original(row_annovar):
    """Varre as colunas de 'Otherinfo' do arquivo do ANNOVAR para obter profundidade e VAF."""
    formato = ""
    valores = ""
    
    colunas_otherinfo = [c for c in row_annovar.index if str(c).startswith('Otherinfo')]
    
    def extrair_numero(nome_col):
        num = re.findall(r'\d+', nome_col)
        return int(num[0]) if num else 0
    colunas_otherinfo.sort(key=extrair_numero)
    
    for idx, col in enumerate(colunas_otherinfo):
        conteudo = str(row_annovar[col]).strip()
        if conteudo.startswith("GT:") or ":GT:" in conteudo or conteudo == "GT":
            formato = conteudo
            if idx + 1 < len(colunas_otherinfo):
                valores = str(row_annovar[colunas_otherinfo[idx + 1]]).strip()
            break
            
    if not formato:
        valores_linha = [str(x).strip() for x in row_annovar.values]
        for i, celula in enumerate(valores_linha):
            if celula.startswith("GT:") or celula == "GT" or ":GT:" in celula:
                formato = celula
                if i + 1 < len(valores_linha):
                    valores = valores_linha[i + 1]
                break

    if not formato or not valores or valores in ['.', 'nan']:
        return '.', '.'
        
    chaves = formato.split(':')
    dados = valores.split(':')
    mapa = {chave: dados[i] for i, chave in enumerate(chaves) if i < len(dados)}
    
    dp = '.'
    vaf = '.'
    
    if 'DP' in mapa and mapa['DP'] != '.':
        dp = mapa['DP']
    elif 'AD' in mapa and mapa['AD'] != '.':
        try:
            reads = [int(x) for x in mapa['AD'].split(',')]
            dp = str(sum(reads))
        except:
            pass
            
    coluna_af = 'AF' if 'AF' in mapa else ('VF' if 'VF' in mapa else None)
    if coluna_af and mapa[coluna_af] != '.':
        try:
            vaf_crua = mapa[coluna_af].split(',')[0]
            val_af = float(vaf_crua)
            if val_af <= 1.0:
                vaf = f"{val_af * 100:.2f}%"
            else:
                vaf = f"{val_af:.2f}%"
        except:
            vaf = mapa[coluna_af]
    elif 'AD' in mapa and mapa['AD'] != '.':
        try:
            ad_valores = mapa['AD'].split(',')
            if len(ad_valores) >= 2:
                ref_count = float(ad_valores[0])
                alt_count = float(ad_valores[1])
                total = ref_count + alt_count
                if total > 0:
                    vaf = f"{(alt_count / total) * 100:.2f}%"
        except:
            pass
            
    return dp, vaf

def extrair_subpopulacao(string_pops, chave):
    if not string_pops or string_pops == '.':
        return '.'
    match = re.search(fr'{chave}:([0-9.]+)', str(string_pops))
    return match.group(1) if match else '.'

def validar_distancia_intronica(hgvs_c):
    if hgvs_c == '.':
        return False
    match = re.search(r'c\.\d+([+-]\d+)', hgvs_c)
    if match:
        try:
            distancia = int(match.group(1))
            if 1 <= distancia <= 5 or -5 <= distancia <= -1:
                return True
        except ValueError:
            pass
    return False

def converter_cancervar(input_cancervar, output_dir, id_amostra):
    print(f"\n📂 Processando arquivo do CancerVar: {input_cancervar}")
    
    # AJUSTE CORRETO: Remove a extensão '.cancervar' para encontrar o ANNOVAR gerado no orquestrador
    input_annovar = input_cancervar.replace('.txt.cancervar', '.txt')
    if not os.path.exists(input_annovar) and input_cancervar.endswith('.cancervar'):
        input_annovar = input_cancervar[:-10] # Remove os 10 caracteres correspondentes a '.cancervar'
        
    df_annovar = None
    if os.path.exists(input_annovar):
        print(f"🔗 [SUCESSO] Arquivo ANNOVAR localizado para cruzamento: {input_annovar}")
        try:
            df_annovar = pd.read_csv(input_annovar, sep='\t', low_memory=False)
            df_annovar.columns = [c.lstrip('#').strip() for c in df_annovar.columns]
        except Exception as e:
            print(f"⚠️ Erro ao ler arquivo do ANNOVAR: {e}")
    else:
        print(f"❌ [ERRO DE CAMINHO] Não achei o arquivo complementar do ANNOVAR em: {input_annovar}")

    try:
        df = pd.read_csv(input_cancervar, sep='\t', low_memory=False)
    except Exception as e:
        print(f"Erro ao ler CancerVar: {e}")
        return

    df.columns = [c.lstrip('#').strip() for c in df.columns]
    
    # ==============================================================================
    # MAPEAMENTO PRÉVIO DE MNVs PARA DEDUPLICAÇÃO
    # ==============================================================================
    posicoes_a_ignorar = set()
    for _, row in df.iterrows():
        chr_val = str(row.get('Chr', '')).strip()
        start_val = str(row.get('Start', '')).strip()
        ref_val = str(row.get('Ref', '')).strip()
        alt_val = str(row.get('Alt', '')).strip()
        
        if not chr_val or chr_val in ['.', '', 'nan'] or 'Chr' in chr_val:
            continue

        if len(ref_val) > 1 and len(ref_val) == len(alt_val) and '-' not in ref_val and '-' not in alt_val:
            try:
                start_num = int(float(start_val))
                for i in range(len(ref_val)):
                    posicoes_a_ignorar.add((chr_val, str(start_num + i)))
            except ValueError:
                pass

    linhas_processadas = []
    cont_sucesso_vaf = 0

    for _, row in df.iterrows():
        chr_val = str(row.get('Chr', '')).strip()
        start_val = str(row.get('Start', '')).strip()
        ref_val = str(row.get('Ref', '')).strip()
        alt_val = str(row.get('Alt', '')).strip()
        
        if not chr_val or chr_val in ['.', '', 'nan'] or 'Chr' in chr_val:
            continue
            
        # FILTRO ATIVO PARA ELIMINAR AS SNVs FRAGMENTADAS REDUNDANTES
        if len(ref_val) == 1 and (chr_val, start_val) in posicoes_a_ignorar:
            continue
            
        func_ref = str(row.get('Func.refGene', '')).strip()
        exonic_func = str(row.get('ExonicFunc.refGene', '')).strip()
        
        # Filtros Biológicos do laboratório
        if 'UTR3' in func_ref or 'UTR5' in func_ref:
            continue
        if 'ncRNA_intronic' in func_ref or 'ncRNA_exonic' in func_ref:
            continue
        if exonic_func == 'synonymous SNV':
            continue

        c_dot, p_dot = extrair_hgvs(row, 'AAChange.refGene')
        
        if func_ref == 'intronic' and not validar_distancia_intronica(c_dot):
            continue

        aachange_bruto = str(row.get('AAChange.refGene', '.'))
        aachange_canonico = aachange_bruto.split(',')[0] if aachange_bruto != '.' else '.'

        # ==============================================================================
        # AJUSTE ROBUSTO INTERVALAR: HERANÇA POR COORDENADA START (MNVs INTEGRADAS)
        # ==============================================================================
        depth, vaf = '.', '.'
        if df_annovar is not None:
            def normalizar_valor(v):
                s = str(v).strip().lower().replace('chr', '')
                if s.endswith('.0'):
                    s = s[:-2]
                return s

            c_search = normalizar_valor(chr_val)
            ref_search = normalizar_valor(ref_val).upper()
            alt_search = normalizar_valor(alt_val).upper()
            
            col_chr = [c for c in df_annovar.columns if 'chr' in c.lower()][0]
            col_start = [c for c in df_annovar.columns if 'start' in c.lower()][0]
            col_ref = [c for c in df_annovar.columns if 'ref' in c.lower() and 'gene' not in c.lower()][0]
            col_alt = [c for c in df_annovar.columns if 'alt' in c.lower()][0]
            
            # Etapa 1: Tentativa de Casamento Perfeito
            match_row = df_annovar[
                (df_annovar[col_chr].apply(normalizar_valor) == c_search) &
                (df_annovar[col_start].apply(normalizar_valor) == normalizar_valor(start_val)) &
                (df_annovar[col_ref].apply(normalizar_valor).str.upper() == ref_search) &
                (df_annovar[col_alt].apply(normalizar_valor).str.upper() == alt_search)
            ]
            
            # Etapa 2: Se falhar (MNVs unificadas), herda os dados da primeira base (coordenada Start)
            if match_row.empty:
                match_posicao_start = df_annovar[
                    (df_annovar[col_chr].apply(normalizar_valor) == c_search) &
                    (df_annovar[col_start].apply(normalizar_valor) == normalizar_valor(start_val))
                ]
                if not match_posicao_start.empty:
                    depth, vaf = extrair_vaf_e_depth_original(match_posicao_start.iloc[0])
            else:
                depth, vaf = extrair_vaf_e_depth_original(match_row.iloc[0])

        if depth != '.' and vaf != '.':
            cont_sucesso_vaf += 1
        
        # Limpezas e padronizações clínicas
        clinvar_raw = str(row.get('clinvar: Clinvar', row.get('clinvar', '.')))
        clinvar_limpo = clinvar_raw.replace('clinvar:', '').strip()
        if clinvar_limpo in ['UNK', '', 'nan', 'None']: clinvar_limpo = '.'

        classificacao_raw = str(row.get('Interp_Prediction', row.get('CancerVar: CancerVar and Evidence', 'Não Classificado')))
        classificacao_somatica_limpa = re.sub(r'^\d+#', '', classificacao_raw).split('EVS=')[0].replace('_', ' ').strip()
        if classificacao_somatica_limpa in ['.', 'nan', '']: classificacao_somatica_limpa = 'Não Classificado'

        pops_gnomad = row.get('Freq_gnomAD_genome_POPs', '.')

        dados_variante = {
            'Chr': chr_val, 'Start': start_val, 'End': row.get('End', '.'), 'Ref': ref_val, 'Alt': alt_val,
            'Gene.refGene': row.get('Ref.Gene', row.get('Gene.refGene', '.')), 'Func.refGene': func_ref,
            'ExonicFunc.refGene': exonic_func, 'AAChange.refGene': aachange_canonico, 'HGVS_cDNA': c_dot,
            'HGVS_Protein': p_dot, 'Depth_Coverage': depth, 'VAF': vaf, 'CancerVar_Classification': classificacao_somatica_limpa,
            'ClinVar_Classification': clinvar_limpo, 
            'esp6500siv2_all': row.get('esp6500siv2_all', '.'), '1000g2015aug_all': row.get('1000g2015aug_all', '.'),
            'ExAC_ALL': row.get('Freq_ExAC_ALL', '.'), 'ExAC_AFR': extrair_subpopulacao(pops_gnomad, 'AFR'),
            'ExAC_AMR': extrair_subpopulacao(pops_gnomad, 'AMR'), 'ExAC_EAS': extrair_subpopulacao(pops_gnomad, 'EAS'),
            'ExAC_FIN': extrair_subpopulacao(pops_gnomad, 'FIN'), 'ExAC_NFE': extrair_subpopulacao(pops_gnomad, 'NFE'),
            'ExAC_OTH': extrair_subpopulacao(pops_gnomad, 'OTH'), 'ExAC_SAS': row.get('ExAC_SAS', '.'),
            'avsnp147': row.get('avsnp147', '.'), 'SIFT_score': row.get('SIFT_score', '.'), 'SIFT_pred': row.get('SIFT_pred', '.'),
            'Polyphen2_HDIV_score': row.get('Polyphen2_HDIV_score', '.'), 'Polyphen2_HDIV_pred': row.get('Polyphen2_HDIV_pred', '.'),
            'Polyphen2_HVAR_score': row.get('Polyphen2_HVAR_score', '.'), 'Polyphen2_HVAR_pred': row.get('Polyphen2_HVAR_pred', '.'),
            'LRT_score': row.get('LRT_score', '.'), 'LRT_pred': row.get('LRT_pred', '.'), 'MutationTaster_score': row.get('MutationTaster_score', '.'),
            'MutationTaster_pred': row.get('MutationTaster_pred', '.'), 'MutationAssessor_score': row.get('MutationAssessor_score', '.'),
            'MutationAssessor_pred': row.get('MutationAssessor_pred', '.'), 'FATHMM_score': row.get('FATHMM_score', '.'),
            'FATHMM_pred': row.get('FATHMM_pred', '.'), 'PROVEAN_score': row.get('PROVEAN_score', '.'), 'PROVEAN_pred': row.get('PROVEAN_pred', '.'),
            'VEST3_score': row.get('VEST3_score', '.'), 'CADD_raw': row.get('CADD_raw', '.'), 'CADD_phred': row.get('CADD_phred', '.'),
            'DANN_score': row.get('DANN_score', '.'), 'fathmm-MKL_coding_score': row.get('fathmm-MKL_coding_score', '.'),
            'fathmm-MKL_coding_pred': row.get('fathmm-MKL_coding_pred', '.'), 'MetaSVM_score': row.get('MetaSVM_score', '.'),
            'MetaSVM_pred': row.get('MetaSVM_pred', '.'), 'MetaLR_score': row.get('MetaLR_score', '.'), 'MetaLR_pred': row.get('MetaLR_pred', '.'),
            'integrated_fitCons_score': row.get('integrated_fitCons_score', '.'), 'integrated_confidence_value': row.get('integrated_confidence_value', '.'),
            'GERP++_RS': row.get('GERP++_RS', '.'), 'phyloP7way_vertebrate': row.get('phyloP7way_vertebrate', '.'),
            'phyloP20way_mammalian': row.get('phyloP20way_mammalian', '.'), 'phastCons7way_vertebrate': row.get('phastCons7way_vertebrate', '.'),
            'phastCons20way_mammalian': row.get('phastCons20way_mammalian', '.'), 'SiPhy_29way_logOdds': row.get('SiPhy_29way_logOdds', '.'),
            'dbscSNV_ADA_SCORE': row.get('dbscSNV_ADA_SCORE', '.'), 'dbscSNV_RF_SCORE': row.get('dbscSNV_RF_SCORE', '.'),
            'dbnsfp31a_interpro': row.get('Interpro_domain', '.'), 'CLNALLELEID': row.get('CLNALLELEID', '.'), 'CLNDN': row.get('CLNDN', '.'),
            'CLNDISDB': row.get('CLNDISDB', '.'), 'CLNREVSTAT': row.get('CLNREVSTAT', '.'), 'CLNSIG': row.get('CLNSIG', '.'),
            'cosmic104': row.get('cosmic91', '.'), 'icgc28': row.get('icgc28', '.'), 'gnomAD_genome_ALL': row.get('Freq_gnomAD_genome_ALL', '.'),
            'gnomAD_genome_AFR': extrair_subpopulacao(pops_gnomad, 'AFR'), 'gnomAD_genome_AMR': extrair_subpopulacao(pops_gnomad, 'AMR'),
            'gnomAD_genome_ASJ': extrair_subpopulacao(pops_gnomad, 'ASJ'), 'gnomAD_genome_EAS': extrair_subpopulacao(pops_gnomad, 'EAS'),
            'gnomAD_genome_FIN': extrair_subpopulacao(pops_gnomad, 'FIN'), 'gnomAD_genome_NFE': extrair_subpopulacao(pops_gnomad, 'NFE'),
            'gnomAD_genome_OTH': extrair_subpopulacao(pops_gnomad, 'OTH')
        }
        linhas_processadas.append(dados_variante)

    print(f"📊 Total de variantes recuperadas com Depth/VAF: {cont_sucesso_vaf} de {len(linhas_processadas)}")

    if not linhas_processadas:
        print("⚠️ Alerta: Nenhuma variante atendeu aos critérios estritos de filtragem biológica.")
        return

    df_final = pd.DataFrame(linhas_processadas)
    colunas_ordenadas = list(dados_variante.keys())
    df_reorganizado = df_final[colunas_ordenadas]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    csv_out = os.path.join(output_dir, f"{id_amostra}_relatorio_final.csv")
    xlsx_out = os.path.join(output_dir, f"{id_amostra}_relatorio_final.xlsx")

    df_reorganizado.to_csv(csv_out, index=False, sep= ';')
    with pd.ExcelWriter(xlsx_out, engine='openpyxl') as writer:
        df_reorganizado.to_excel(writer, index=False, sheet_name='Variantes_Somáticas')
        
    print(f"✅ Relatório criado com sucesso em: {xlsx_out}")

if __name__ == "__main__":
    if len(sys.argv) == 4:
        converter_cancervar(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        # Configuração para execução local/manual direta na pasta de anotação
        pasta_anotacao = "/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/05_anotacao"
        padrao_busca = os.path.join(pasta_anotacao, "*_cancervar.output.hg38_multianno.txt.cancervar")
        arquivos_encontrados = glob.glob(padrao_busca)
        
        if arquivos_encontrados:
            input_padrao = arquivos_encontrados[0]
            # Coleta o nome base da amostra de forma correta (ex: 1204)
            id_amostra_padrao = os.path.basename(input_padrao).split('_cancervar')[0]
            output_dir_padrao = "/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/06_relatorios_finais"
            converter_cancervar(input_padrao, output_dir_padrao, id_amostra_padrao)
        else:
            print("❌ Nenhum arquivo .cancervar correspondente ao padrão foi localizado.")
