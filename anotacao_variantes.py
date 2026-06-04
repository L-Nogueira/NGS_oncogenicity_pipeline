import os
import subprocess

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
    
    prefixo_saida = os.path.join(pasta_saida_abs, f"{id_amostra}_annotated")
    tabela_annovar_txt = f"{prefixo_saida}.hg38_multianno.txt"

    print(f"🧬 [Anotação] Iniciando processo para a amostra: {id_amostra}")

    # ==========================================
    # PASSO 1: Anotação ANNOVAR Local
    # ==========================================
    cmd_annovar = (
        f"perl {TABLE_ANNOVAR} {os.path.abspath(vcf_entrada)} {HUMANDB} "
        f"-buildver hg38 "
        f"-outfile {prefixo_saida} "
        f"-remove "
        f"-protocol refGene,clinvar_20240917 "
        f"-operation g,f "
        f"-nastring . "
        f"-otherinfo"
    )

    try:
        print(f"🏃‍♂️ [Passo 1/2] Executando ANNOVAR (Gene + ClinVar)...")
        subprocess.run(cmd_annovar, shell=True, check=True)
        print(f"✅ ANNOVAR concluído com sucesso.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Erro na anotação do ANNOVAR: {e}")
        return False

    # ==========================================
    # PASSO 2: Classificação Dinâmica (CancerVar)
    # ==========================================
    tabela_annovar_txt_abs = os.path.abspath(tabela_annovar_txt)
    arquivo_cancervar_saida_abs = os.path.join(pasta_saida_abs, f"{id_amostra}_cancervar.output")
    cancervar_config_abs = os.path.abspath(cancervar_config)
    cancervar_py_abs = os.path.abspath(cancervar_py)
    
    # Descobre o diretório raiz do CancerVar para isolar o ambiente de execução
    dir_cancervar_raiz = os.path.dirname(cancervar_py_abs)
    cancervardb_abs = os.path.join(dir_cancervar_raiz, "cancervardb")

    # Montagem limpa do comando
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
        print(f"🏃‍♂️ [Passo 2/2] Executando predição dinâmica do CancerVar...")
        
        # A MUDANÇA ESSENCIAL: cwd=dir_cancervar_raiz força o script a rodar de dentro da 
        # pasta do próprio CancerVar. Isso faz com que as checagens internas de pastas relativas 
        # funcionem sem precisar de links simbólicos locais no seu diretório de projetos.
        subprocess.run(cmd_cancervar, shell=True, check=True, cwd=dir_cancervar_raiz)
        
        print(f"✅ CancerVar concluído com sucesso!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Erro na execução do CancerVar: {e}")
        return False
