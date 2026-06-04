import os
import glob
import subprocess
import sys
import time
import re

# ======================================================================
# --- CONFIGURAÇÃO DE CAMINHOS E EXECUTÁVEIS ---
# ======================================================================
PASTA_BRUTOS = "dados_brutos"
REFERENCIA = os.path.expanduser("~/laboratorio_bioinfo/genomas_referencia/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
PASTA_VARIANTES = "04_variantes"
PASTA_ANOTACAO = "05_anotacao" 

# Arquivos e caminhos de suporte
ARQUIVO_BED = "/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/arquivo .bed/targeted_regions_224441_ensembl.bed"
ANNOVAR_DIR = "/home/l.nogueira/laboratorio_bioinfo/softwares/annovar"

# CORREÇÃO CANCERVAR: Centralizando caminhos e corrigindo flags de argumento
CANCERVAR_DIR = "/home/l.nogueira/laboratorio_bioinfo/softwares/CancerVar"
CANCERVAR_PY = os.path.join(CANCERVAR_DIR, "CancerVar.py")
CANCERVAR_CONFIG = os.path.join(CANCERVAR_DIR, "config.ini")

# Caminho do seu novo script de conversão amigável
SCRIPT_CONVERSOR = "/home/l.nogueira/laboratorio_bioinfo/scripts_bioinfo/converter_output.py"

# ======================================================================

def executar_script(comando):
    """Executa scripts simples da pipeline (triagem, limpeza, variantes)."""
    try:
        subprocess.run(comando, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Erro crítico no passo: {comando}\nInterrompendo a esteira.")
        sys.exit(1)

def contar_reads_fastq(caminho_fastq):
    """Conta rapidamente o número aproximado de reads dividindo o total de linhas por 4."""
    try:
        cmd = f"zcat '{caminho_fastq}' | wc -l"
        linhas = int(subprocess.check_output(cmd, shell=True).decode().strip())
        return linhas // 4
    except:
        return 5000000 # Valor padrão estimado caso a contagem falhe

def exibir_barra(atual, total, prefixo=''):
    """Gera uma barra de progresso visual no terminal."""
    largura = 40
    progresso = int(largura * atual / total)
    barra = '█' * progresso + '-' * (largura - progresso)
    percentual = round((atual / total) * 100, 1)
    print(f"\r{prefixo} |{barra}| {percentual}% ({atual}/{total})", end='', flush=True)

def rodar_alinhamento_com_progresso(r1, r2, referencia, id_amostra, total_reads_esperado):
    """Executa o BWA-MEM e Samtools monitorando o arquivo SAM temporário para exibir a barra."""
    import alinhamento
    alinhamento.realizar_alinhamento(r1, r2, referencia)

# --- INÍCIO DA ESTEIRA ---
print("======================================================================")
print("🧬  ESTEIRA DE BIOINFORMÁTICA NGS AUTOMATIZADA COM PONTOS DE CHECAGEM  🧬")
print("======================================================================\n")

tempo_inicio_global = time.time()

# Buscar pares Forward (R1) na pasta de dados brutos
arquivos_r1 = sorted(glob.glob(os.path.join(PASTA_BRUTOS, "*_R1_*.fastq.gz")))

if not arquivos_r1:
    print(f"❌ Nenhum arquivo R1 encontrado na pasta '{PASTA_BRUTOS}'. Verifique o caminho.")
    sys.exit(1)

total_amostras = len(arquivos_r1)
print(f"📋 Foram encontradas {total_amostras} amostras para processar.\n")

for index, r1_bruto in enumerate(arquivos_r1, start=1):
    tempo_inicio_amostra = time.time()
    
    # Deduzir o arquivo R2 correspondente
    r2_bruto = r1_bruto.replace("_R1_", "_R2_")
    if not os.path.exists(r2_bruto):
        print(f"⚠️ Alerta: Par R2 não encontrado para {r1_bruto}. Pulando amostra.")
        continue

    # Identificar o ID da amostra
    nome_base_r1 = os.path.basename(r1_bruto)
    id_amostra = nome_base_r1.split('_R1')[0]

    print("\n" + "="*70)
    print(f"🔬 PROCESSANDO AMOSTRA ({index}/{total_amostras}): {id_amostra}")
    print("="*70)

    # --- PASSO 1: TRIAGEM DE QUALIDADE (FASTQC) ---
    fastqc_html_r1 = os.path.join("fastqc_results", nome_base_r1.replace(".fastq.gz", "_fastqc.html"))
    fastqc_html_r2 = os.path.join("fastqc_results", os.path.basename(r2_bruto).replace(".fastq.gz", "_fastqc.html"))

    if os.path.exists(fastqc_html_r1) and os.path.exists(fastqc_html_r2):
        print("⏭️  [PASSO 1] FastQC já executado anteriormente. Pulando...")
    else:
        print("📊 [PASSO 1] Executando triagem de integridade e FastQC...")
        executar_script(f"python3 ../scripts_bioinfo/triagem.py {r1_bruto} {r2_bruto}")

    # --- PASSO 2: FILTRAGEM E LIMPEZA (TRIMMOMATIC) ---
    r1_trimado = os.path.join("trimmed_data", nome_base_r1.replace(".fastq.gz", "_trim_P.fastq.gz"))
    r2_trimado = os.path.join("trimmed_data", os.path.basename(r2_bruto).replace(".fastq.gz", "_trim_P.fastq.gz"))

    if os.path.exists(r1_trimado) and os.path.exists(r2_trimado):
        print("⏭️  [PASSO 2] Dados já filtrados pelo Trimmomatic. Pulando...")
    else:
        print("✂️  [PASSO 2] Executando limpeza de adaptadores e qualidade baixa...")
        executar_script(f"python3 ../scripts_bioinfo/limpeza.py {r1_bruto} {r2_bruto}")

    # --- PASSO 3: ALINHAMENTO CONTRA O GENOMA (BWA-MEM) ---
    id_alinhamento = id_amostra
    bam_ordenado = os.path.join("alignment_results", f"{id_alinhamento}_sorted.bam")

    if os.path.exists(bam_ordenado):
        print("⏭️  [PASSO 3] Alinhamento BAM já existe e está ordenado. Pulando...")
    else:
        print("🚀 [PASSO 3] Calculando tamanho da amostra para mapeamento...")
        total_reads_amostra = contar_reads_fastq(r1_trimado)
        print(f"📉 Total de leituras estimadas para processamento: {total_reads_amostra}")
        print("🧬 Mapeando reads contra o genoma de referência...")
        rodar_alinhamento_com_progresso(r1_trimado, r2_trimado, REFERENCIA, id_alinhamento, total_reads_amostra)

    # --- PASSO 4: CHAMADA DE VARIANTES SOMÁTICAS (MUTECT2) ---
    vcf_final = os.path.join(PASTA_VARIANTES, f"{id_amostra}_variants.vcf")
    
    if os.path.exists(vcf_final):
        print(f"⏭️  [PASSO 4] Variantes somáticas (VCF) já calculadas para {id_amostra}. Pulando...")
    else:
        print(f"📦 [PASSO 4] Iniciando Mutect2 Smart Target para a amostra {id_amostra}...")
        
        sys.path.append('../scripts_bioinfo')
        import variantes_mutect2

        sucesso = variantes_mutect2.rodar_mutect2_smart_target(
            id_amostra=id_amostra,
            referencia=REFERENCIA,
            bam_entrada=bam_ordenado,
            arquivo_bed=ARQUIVO_BED,
            pasta_saida=PASTA_VARIANTES
        )

        if not sucesso:
            print(f"❌ Falha no processamento somático da amostra {id_amostra}")
            continue

    # --- PASSO 5: ANOTAÇÃO (CLINVAR + CANCERVAR) ---
    arquivo_check_anotacao = os.path.join(PASTA_ANOTACAO, f"{id_amostra}_cancervar.output")
    
    if os.path.exists(arquivo_check_anotacao):
        print(f"⏭️  [PASSO 5] Anotação Clínica e Oncogenética já realizada. Pulando...")
    else:
        print(f"🧬 [PASSO 5] Iniciando Anotação (ClinVar + CancerVar) para {id_amostra}...")
        
        sys.path.append('../scripts_bioinfo')
        import anotacao_variantes

        # Altera temporariamente para a pasta do CancerVar para evitar quebra de caminhos internos do software
        diretorio_original = os.getcwd()
        try:
            os.chdir(CANCERVAR_DIR)
            sucesso_anotacao = anotacao_variantes.rodar_anotacao(
                id_amostra=id_amostra,
                vcf_entrada=os.path.join(diretorio_original, vcf_final),
                pasta_saida=os.path.join(diretorio_original, PASTA_ANOTACAO),
                cancervar_py=CANCERVAR_PY,      
                cancervar_config=CANCERVAR_CONFIG 
            )
        finally:
            # Retorna com segurança para a pasta original do projeto
            os.chdir(diretorio_original)

        if not sucesso_anotacao:
            print(f"⚠️ Aviso: A anotação encontrou problemas para a amostra {id_amostra}")
        else:
            print("✅ CancerVar concluído com sucesso!")

    # --- PASSO EXTRA AUTOMATIZADO: CONVERSÃO E ORGANIZAÇÃO SEPARADA POR AMOSTRA ---
    PASTA_RELATORIOS = "06_relatorios_finais"
    cancervar_final_txt = os.path.join(PASTA_ANOTACAO, f"{id_amostra}_cancervar.output.hg38_multianno.txt.cancervar")
    
    # Define a pasta exclusiva desta amostra dentro do diretório final de relatórios
    pasta_destino_amostra = os.path.join(PASTA_RELATORIOS, id_amostra)

    if os.path.exists(cancervar_final_txt):
        print(f"📊 [PASSO EXTRA] Formatando e movendo tabelas de {id_amostra} para pasta de relatórios dedicada...")
        # Executa o script passando: 1. O arquivo txt bruto, 2. A nova pasta de destino, 3. O ID da amostra
        comando_conversao = f"python3 '{SCRIPT_CONVERSOR}' '{cancervar_final_txt}' '{pasta_destino_amostra}' '{id_amostra}'"
        executar_script(comando_conversao)
    else:
        print(f"⚠️ Alerta: Arquivo final .cancervar não localizado para conversão.")
    # Métricas de tempo
    tempo_de_corrida_amostra = time.time() - tempo_inicio_amostra
    tempo_parcial_global = time.time() - tempo_inicio_global

    print("\n" + "-"*70)
    print(f"✅ Amostra {id_amostra} concluída/verificada!")
    print(f"⏱️  Tempo desta amostra: {tempo_de_corrida_amostra/60:.2f} minutos")
    print(f"⏱️  Tempo total acumulado da esteira: {tempo_parcial_global/60:.2f} minutos")
    print("-"*70 + "\n")

tempo_total_global = time.time() - tempo_inicio_global
print("======================================================================")
print(f"🎉 ESTEIRA CONCLUÍDA COM SUCESSO EM {tempo_total_global/60:.2f} MINUTOS!")
print("======================================================================")
