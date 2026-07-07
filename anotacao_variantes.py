import os
import subprocess
import pandas as pd
import re

def pre_filtrar_tabela_annovar(caminho_txt):
    """
    Filtra o arquivo gerado pelo ANNOVAR mantendo apenas variantes exônicas 
    e de splicing de interesse, descartando UTRs, íntrons profundos e sinônimas.
    """
    print("🧹 [Pré-Filtro] Limpando variantes intrônicas e UTRs para acelerar o CancerVar...")
    
    # Lê a tabela do ANNOVAR
    df = pd.read_csv(caminho_txt, sep='\t', low_memory=False)
    df.columns = [c.lstrip('#').strip() for c in df.columns]
    
    def validar_distancia_intronica(aachange):
        val = str(aachange)
        if not val or val == '.' or val == 'nan': 
            return False
        primeiro = val.split(',')[0]
        match = re.search(r'c\.\d+([+-]\d+)', primeiro)
        if match:
            try:
                dist = int(match.group(1))
                return 1 <= dist <= 5 or -5 <= dist <= -1
            except ValueError: 
                pass
        return False

    linhas_filtradas = []
    for idx, row in df.iterrows():
        func_ref = str(row.get('Func.refGene', '')).strip()
        exonic_func = str(row.get('ExonicFunc.refGene', '')).strip()
        aa_change = str(row.get('AAChange.refGene', '')).strip()
        
        # Filtros biológicos idênticos aos do conversor final
        if 'UTR3' in func_ref or 'UTR5' in func_ref: 
            continue
        if 'ncRNA_intronic' in func_ref or 'ncRNA_exonic' in func_ref: 
            continue
        if exonic_func == 'synonymous SNV': 
            continue
        if func_ref == 'intronic' and not validar_distancia_intronica(aa_change): 
            continue
        
        linhas_filtradas.append(row)
        
    # Recria o arquivo .txt substituindo pelo conteúdo filtrado
    df_filtrado = pd.DataFrame(linhas_filtradas)
    df_filtrado.to_csv(caminho_txt, sep='\t', index=False)
    print(f"降低 Tabela enxugada de {len(df)} para {len(df_filtrado)} variantes de interesse somático.")

def rodar_anotacao(id_amostra, vcf_entrada, pasta_saida, cancervar_py, cancervar_config):
    """
    Realiza a anotação de variantes utilizando ANNOVAR (com ClinVar atualizado)
    e complementa com a classificação dinâmica do CancerVar usando a tabela gerada.
    """
    # Garantir caminhos absolutos robustos para as pastas do projeto
    pasta_saida_abs = os.path.abspath(pasta_saida)
    os.makedirs(pasta_saida_abs, exist_ok=True)
    
    ANNOVAR_DIR = os.path.expanduser("~/laboratorio_bioinfo/softwares/annovar")
    HUMANDB = os.path.join(ANNOVAR_DIR, "humandb")
    TABLE_ANNOVAR = os.path.join(ANNOVAR_DIR, "table_annovar.pl")
    
    # CORREÇÃO: Nome limpo sem o sufixo '_franklin' para compatibilidade nativa com o CancerVar
    prefixo_saida = os.path.join(pasta_saida_abs, id_amostra)
    tabela_annovar_txt = f"{prefixo_saida}.hg38_multianno.txt"
    vcf_annovar_gerado = f"{prefixo_saida}.hg38_multianno.vcf"

    print(f"🧬 [Anotação] Iniciando processo para a amostra: {id_amostra}")

    # ==========================================
    # PASSO 1: ANNOVAR Adaptado (Gera a estrutura enxuta em segundos)
    # ==========================================
    # Ajustado o output para o padrão que o CancerVar espera receber e pular buscas pesadas
    cmd_annovar = (
        f"perl {TABLE_ANNOVAR} {os.path.abspath(vcf_entrada)} {HUMANDB} "
        f"-buildver hg38 -outfile {prefixo_saida} -remove "
        f"-protocol refGene,clinvar_20240917,cosmic91 -operation g,f,f "
        f"-nastring . -vcfinput"
    )

    try:
        # Checagem inteligente: se o arquivo final do CancerVar já existe, 
        # o ANNOVAR inicial também pode ser pulado com segurança total.
        arquivo_cancervar_saida_abs = os.path.join(pasta_saida_abs, f"{id_amostra}_cancervar.output")
        cancervar_final_txt = arquivo_cancervar_saida_abs + ".hg38_multianno.txt.cancervar"

        if not os.path.exists(cancervar_final_txt):
            if not os.path.exists(tabela_annovar_txt):
                print(f"🏃‍♂️ [Passo 1/2] Executando ANNOVAR Inicial (TXT + VCF)...")
                subprocess.run(cmd_annovar, shell=True, check=True)
                print("✅ Etapa ANNOVAR inicial concluída.")
            else:
                print("⏭️  [Passo 1/2] Arquivo base do ANNOVAR já localizado.")
        else:
            print("⏭️  [SKIP GLOBAL] Todos os outputs de Anotação/CancerVar já existem para esta amostra.")
            return True

    except Exception as e:
        print(f"❌ Erro crítico ao executar o ANNOVAR: {e}")
        return False

    # ==========================================
    # PASSO INTERMEDIÁRIO: Filtragem da Tabela ANNOVAR
    # ==========================================
    if os.path.exists(tabela_annovar_txt):
        pre_filtrar_tabela_annovar(tabela_annovar_txt)

    # ==========================================
    # PASSO 2: Predição Dinâmica do CancerVar
    # ==========================================
    tabela_annovar_txt_abs = os.path.abspath(tabela_annovar_txt)
    arquivo_cancervar_saida_abs = os.path.join(pasta_saida_abs, f"{id_amostra}_cancervar.output")
    cancervar_config_abs = os.path.abspath(cancervar_config)
    cancervar_py_abs = os.path.abspath(cancervar_py)
    
    dir_cancervar_raiz = os.path.dirname(cancervar_py_abs)
    cancervardb_abs = os.path.join(dir_cancervar_raiz, "cancervardb")

    cmd_cancervar = (
        f"python3 {cancervar_py_abs} "
        f"-c {cancervar_config_abs} "
        f"-i {tabela_annovar_txt_abs} "
        f"-t AVinput "
        f"-d {cancervardb_abs} "
        f"-b hg38 "
        f"-o {arquivo_cancervar_saida_abs}"
    )

    try:
        # Checar se o output final do CancerVar já existe para permitir o SKIP inteligente por amostra
        cancervar_final_txt = arquivo_cancervar_saida_abs + ".hg38_multianno.txt.cancervar"
        if not os.path.exists(cancervar_final_txt):
            print(f"🏃‍♂️ [Passo 2/2] Executando predição dinâmica do CancerVar...")
            # Mudança de diretório de trabalho (cwd) necessária devido às dependências internas do CancerVar
            subprocess.run(cmd_cancervar, shell=True, check=True, cwd=dir_cancervar_raiz)
            print("✅ CancerVar concluído com sucesso!")
        else:
            print("⏭️  [Passo 2/2] Output do CancerVar já existe. Pulando...")
        return True
            
    except Exception as e:
        print(f"❌ Erro crítico ao executar o CancerVar: {e}")
        return False
