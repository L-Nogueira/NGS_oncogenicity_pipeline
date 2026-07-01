import pandas as pd
import sys
import os
import re
import glob

def extrair_hgvs(row, coluna_fonte):
    """Função auxiliar para extrair cDNA (c.) e Proteína (p.) da string do ANNOVAR."""
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

def extrair_vaf_e_depth(row):
    """
    Varre cirurgicamente as colunas nominais de 'Otherinfo' geradas pelo ANNOVAR/CancerVar.
    Localiza onde está a chave 'GT:' e pareia com a coluna de valores subsequente.
    """
    formato = ""
    valores = ""
    
    # 1. Identifica quais colunas do DataFrame começam com 'Otherinfo'
    colunas_otherinfo = [c for c in row.index if str(c).startswith('Otherinfo')]
    
    # Ordena as colunas numericamente (Otherinfo1, Otherinfo2... Otherinfo10, Otherinfo11)
    def extrair_numero(nome_col):
        num = re.findall(r'\d+', nome_col)
        return int(num[0]) if num else 0
    colunas_otherinfo.sort(key=extrair_numero)
    
    # 2. Procura pela coluna que dita o formato (contém 'GT:')
    for idx, col in enumerate(colunas_otherinfo):
        conteudo = str(row[col]).strip()
        if conteudo.startswith("GT:") or ":GT:" in conteudo or conteudo == "GT":
            formato = conteudo
            # No padrão do CancerVar, os valores numéricos estão estritamente na coluna Otherinfo seguinte
            if idx + 1 < len(colunas_otherinfo):
                valores = str(row[colunas_otherinfo[idx + 1]]).strip()
            break
            
    # Se falhar no mapeamento básico por nome de coluna, faz o fallback linear por conteúdo da linha
    if not formato:
        valores_linha = [str(x).strip() for x in row.values]
        for i, celula in enumerate(valores_linha):
            if celula.startswith("GT:") or celula == "GT":
                formato = celula
                if i + 1 < len(valores_linha):
                    valores = valores_linha[i + 1]
                break

    # Se não localizou a estrutura clássica de VCF, retorna ponto de segurança
    if not formato or not valores or valores == '.' or valores == 'nan':
        return '.', '.'
        
    chaves = formato.split(':')
    dados = valores.split(':')
    
    # Cria o dicionário pareando chave-valor por posição
    mapa = {chave: dados[i] for i, chave in enumerate(chaves) if i < len(dados)}
    
    dp = '.'
    vaf = '.'
    
    # 3. Extração precisa da Cobertura Total (DP)
    if 'DP' in mapa and mapa['DP'] != '.':
        dp = mapa['DP']
    elif 'AD' in mapa and mapa['AD'] != '.':
        try:
            # Fallback caso o DP não esteja explícito: soma reads de Ref + Alt
            reads = [int(x) for x in mapa['AD'].split(',')]
            dp = str(sum(reads))
        except:
            pass
            
    # 4. Extração precisa da Frequência Alélica / VAF (AF ou VF)
    coluna_af = 'AF' if 'AF' in mapa else ('VF' if 'VF' in mapa else None)
    if coluna_af and mapa[coluna_af] != '.':
        try:
            # Pega o primeiro valor (tratando possíveis multialélicos separados por vírgula)
            vaf_crua = mapa[coluna_af].split(',')[0]
            val_af = float(vaf_crua)
            
            # Correção matemática para interpretar frequências como 1.0000 ou 0.2500 corretamente
            if val_af <= 1.0:
                vaf = f"{val_af * 100:.2f}%"
            else:
                # Caso o variant caller já envie em escala de 0-100 (ex: 25.4)
                vaf = f"{val_af:.2f}%"
        except:
            vaf = mapa[coluna_af]
            
    return dp, vaf

def filtrar_regras_laboratorio(row):
    """Filtro Clínico Restrito do Laboratório (Apenas Exônicos válidos e Íntrons 1-5)."""
    func = str(row.get('Func.refGene', '')).strip()
    aachange = str(row.get('AAChange.refGene', '')).strip()
    
    if func == 'intronic':
        detalhe = str(row.get('GeneDetail.refGene', '')).strip()
        if not detalhe or detalhe == '.' or detalhe == 'nan':
            return False
        matches = re.findall(r'[c|*]\.[0-9]+([+-])([0-9]+)', detalhe)
        if not matches:
            return False
        for sinal, numero_str in matches:
            distancia = int(numero_str)
            if 1 <= distancia <= 5:
                return True
        return False

    if not aachange or aachange == '.' or aachange == 'nan':
        return False
        
    return True

def converter_cancervar(input_txt, output_dir, id_amostra):
    print(f"📊 Lendo arquivo bruto do CancerVar: {os.path.basename(input_txt)}")
    
    try:
        df = pd.read_csv(input_txt, sep='\t', low_memory=False)
        
        if df.iloc[0, 1] == 'Start':
            df = df.drop(df.index[0]).reset_index(drop=True)
            
        df.columns = [c.strip() for c in df.columns]
        
        # 1. Aplicação do Filtro do Laboratório (Remove ncRNAs e UTRs sem impacto)
        df = df[df.apply(filtrar_regras_laboratorio, axis=1)].reset_index(drop=True)

        if len(df) == 0:
            print("⚠️ Aviso: Nenhuma variante restou após os filtros de região.")
            return

        # 2. Ajuste de Múltiplos Transcritos
        for col_aachange in ['AAChange.refGene', 'AAChange.ensGene', 'AAChange.knownGene']:
            if col_aachange in df.columns:
                df[col_aachange] = df[col_aachange].astype(str).apply(
                    lambda x: x.split(',')[0] if x and x != '.' and x != 'nan' else x
                )
        
        if df.columns[0].startswith('CancerVar:'):
            df.rename(columns={df.columns[0]: 'Interpretation_Details'}, inplace=True)

        # 3. Extração Automática de HGVS
        if 'AAChange.refGene' in df.columns:
            hgvs_extraido = df.apply(lambda row: extrair_hgvs(row, 'AAChange.refGene'), axis=1)
            df['HGVS_cDNA'] = [h[0] for h in hgvs_extraido]
            df['HGVS_Protein'] = [h[1] for h in hgvs_extraido]
        else:
            df['HGVS_cDNA'] = '.'
            df['HGVS_Protein'] = '.'

        def preencher_hgvs_intron(row):
            if row['HGVS_cDNA'] == '.' and row['Func.refGene'] == 'intronic':
                detalhe = str(row.get('GeneDetail.refGene', ''))
                match = re.search(r'(c\.[^,:\s]+)', detalhe)
                if match:
                    return match.group(1)
            return row['HGVS_cDNA']
        df['HGVS_cDNA'] = df.apply(preencher_hgvs_intron, axis=1)

        # ------------------------------------------------------------------
        # AJUSTE EXATO: EXTRAÇÃO VIA PAR CHAVE-VALOR DO CANCERVAR
        # ------------------------------------------------------------------
        print("📈 Extraindo Depth e VAF pareando as colunas Otherinfo...")
        metricas = df.apply(extrair_vaf_e_depth, axis=1)
        df['Depth_Coverage'] = [m[0] for m in metricas]
        df['VAF'] = [m[1] for m in metricas]
        # ------------------------------------------------------------------

        if 'Interpretation_Details' in df.columns:
            df['CancerVar_Classification'] = df['Interpretation_Details'].astype(str).apply(
                lambda x: x.split('Link:')[0].strip() if 'Link:' in x else x.strip()
            )
        else:
            df['CancerVar_Classification'] = 'Não Classificado'

    except Exception as e:
        print(f"❌ Erro ao ler e processar o arquivo: {e}")
        return

    # Garante a remoção completa de todas as colunas "Otherinfo" originais e poluídas
    colunas_finais = [c for c in df.columns if not str(c).startswith('Otherinfo')]

    colunas_para_o_inicio = [
        'Chr', 'Start', 'End', 'Ref', 'Alt', 
        'Gene.refGene', 'Func.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
        'HGVS_cDNA', 'HGVS_Protein', 'Depth_Coverage', 'VAF', 'CancerVar_Classification'
    ]
    
    colunas_para_o_inicio = [c for c in colunas_para_o_inicio if c in colunas_finais]
    outras_colunas = [c for c in colunas_finais if c not in colunas_para_o_inicio]
    df_reorganizado = df[colunas_para_o_inicio + outras_colunas]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    csv_out = os.path.join(output_dir, f"{id_amostra}_relatorio_final.csv")
    xlsx_out = os.path.join(output_dir, f"{id_amostra}_relatorio_final.xlsx")

    df_reorganizado.to_csv(csv_out, index=False, sep=';')
    with pd.ExcelWriter(xlsx_out, engine='openpyxl') as writer:
        df_reorganizado.to_excel(writer, index=False, sheet_name='Variantes_Somáticas')
    print(f"✅ Sucesso absoluto! Nova Planilha Excel gerada com Cobertura e VAF em: {xlsx_out}")

if __name__ == "__main__":
    if len(sys.argv) == 4:
        input_file = sys.argv[1]
        output_dir = sys.argv[2]
        id_amostra = sys.argv[3]
        converter_cancervar(input_file, output_dir, id_amostra)
    else:
        pasta_anotacao = "/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/05_anotacao"
        padrao_busca = os.path.join(pasta_anotacao, "*cancervar.output.hg38_multianno.txt")
        arquivos_encontrados = glob.glob(padrao_busca)
        
        if arquivos_encontrados:
            input_padrao = arquivos_encontrados[0]
            nome_arquivo = os.path.basename(input_padrao)
            id_amostra_padrao = nome_arquivo.split('_cancervar')[0]
            dir_out_padrao = os.path.join("/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/06_relatorios_finais", id_amostra_padrao)
            converter_cancervar(input_padrao, dir_out_padrao, id_amostra_padrao)
        else:
            print(f"⚠️ Nenhum arquivo de anotação localizado.")
