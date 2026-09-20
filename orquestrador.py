import os
import glob
import subprocess
import sys
import time
import re

# Módulos customizados do pipeline
import triagem
import limpeza
import alinhamento
import variantes_mutect2
import anotacao_variantes

# ======================================================================
# --- CONFIGURAÇÃO DE CAMINHOS E EXECUTÁVEIS (PORTABILIDADE) ---
# ======================================================================

# Captura dinamicamente a Home do usuário:
HOME = os.path.expanduser("~")

#Caminhos de diretórios internos
PASTA_BRUTOS = "dados_brutos"
PASTA_VARIANTES = "04_variantes"
PASTA_ANOTACAO = "05_anotacao" 

# Genoma de Referência
REFERENCIA = os.path.join(HOME, "laboratorio_bioinfo", "genomas_referencia", "Homo_sapiens.GRCh38.dna.primary_assembly.fa")

# BED
ARQUIVO_BED = os.path.join(HOME, "laboratorio_bioinfo", "projetos_miseq_real", "arquivo .bed", "targeted_regions_FIXED.bed")

# ANNOVAR / CANCERVAR
ANNOVAR_DIR = os.path.join(HOME, "laboratorio_bioinfo", "softwares", "annovar")
CANCERVAR_DIR = os.path.join(HOME, "laboratorio_bioinfo", "softwares", "CancerVar")
CANCERVAR_PY = os.path.join(CANCERVAR_DIR, "CancerVar.py")
CANCERVAR_CONFIG = os.path.join(CANCERVAR_DIR, "config.ini")

# Caminho do script de conversão (Excel/CSV)
SCRIPT_CONVERSOR = os.path.join(HOME, "laboratorio_bioinfo", "scripts_bioinfo", "converter_output.py")
# ======================================================================

def executar_script(comando):
    """Função auxiliar para executar comandos de terminal monitorando erros."""
    try:
        subprocess.run(comando, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Erro crítico na execução do comando: {comando}\nRetorno: {e}")
        sys.exit(1)

def obter_id_amostra(arquivo_r1):
    """Extrai o ID da amostra isolando o padrão antes do primeiro _R1."""
    base = os.path.basename(arquivo_r1)
    match = re.match(r"^([A-Za-z0-9_-]+)_R1", base)
    if match:
        return match.group(1)
    return base.split('_R1')[0]

# ======================================================================
# --- EXECUÇÃO DO PIPELINE ORQUESTRADO ---
# ======================================================================

if __name__ == "__main__":
    tempo_inicio_global = time.time()
    print("======================================================================")
    print("🧬 INICIANDO PIPELINE AUTOMATIZADO DE BIOINFORMÁTICA (MiSeq Real) 🧬")
    print("======================================================================")

    # Criar pastas essenciais caso não existam
    os.makedirs(PASTA_VARIANTES, exist_ok=True)
    os.makedirs(PASTA_ANOTACAO, exist_ok=True)

    # --- NOVIDADE PARA PARALELIZAÇÃO ---
    # Verifica se o ID da amostra foi passado diretamente como argumento via terminal
    if len(sys.argv) > 1:
        ids_amostras = [sys.argv[1]]
        print(f"🚀 [Modo Paralelo] Executando especificamente a amostra: {ids_amostras[0]}")
    else:
        # Se você rodar sem passar argumentos, ele mantém o comportamento antigo de listar tudo
        print("🧬 [Modo Sequencial] Detectando todas as amostras no sistema...")
        vcf_mutect_lista = glob.glob(os.path.join(PASTA_VARIANTES, "*_variants.vcf"))
        vcf_dragen_lista = glob.glob(os.path.join(PASTA_VARIANTES, "*.hard-filtered.vcf"))
        arquivos_r1 = glob.glob(os.path.join(PASTA_BRUTOS, "*_R1_*.fastq.gz"))
        
        ids_amostras = set()
        
        for vcf_path in vcf_mutect_lista:
            id_vcf = os.path.basename(vcf_path).replace("_variants.vcf", "")
            ids_amostras.add(id_vcf)
            
        for vcf_path in vcf_dragen_lista:
            id_vcf = os.path.basename(vcf_path).replace(".hard-filtered.vcf", "")
            ids_amostras.add(id_vcf)
            
        for r1_path in arquivos_r1:
            ids_amostras.add(obter_id_amostra(r1_path))
            
        ids_amostras = sorted(list(ids_amostras))

    if not ids_amostras:
        print(f"❌ Erro: Nenhuma amostra localizada.")
        sys.exit(1)

    if len(ids_amostras) > 1:
        print(f"📦 Total de amostras identificadas para processamento sequencial: {len(ids_amostras)}")

    # ======================================================================
    # O LOOP ABAIXO SEGUE EXATAMENTE IGUAL AO QUE JÁ ESTAVA FUNCIONANDO
    # ======================================================================
    for id_amostra in ids_amostras:
        tempo_inicio_amostra = time.time()
        
        print("\n" + "-"*80)
        print(f"🔬 PROCESSANDO AMOSTRA: {id_amostra}")
        print(f"{'-'*80}")

        # Definir caminhos possíveis para os arquivos VCF desta amostra
        vcf_mutect2 = os.path.join(PASTA_VARIANTES, f"{id_amostra}_variants.vcf")
        vcf_dragen = os.path.join(PASTA_VARIANTES, f"{id_amostra}.hard-filtered.vcf")
        
        # Atribuição dinâmica do VCF de entrada real baseado no que existir na pasta
        if os.path.exists(vcf_mutect2):
            vcf_entrada_real = vcf_mutect2
        elif os.path.exists(vcf_dragen):
            vcf_entrada_real = vcf_dragen
        else:
            vcf_entrada_real = vcf_mutect2

        # DEFINIÇÃO DE MODO: Se algum dos dois VCFs já existe, pula upstream
        if os.path.exists(vcf_mutect2) or os.path.exists(vcf_dragen):
            print(f"ℹ️  Modo VCF Direto detectado para {id_amostra} ({os.path.basename(vcf_entrada_real)} encontrado).")
            print(f"⏩ [SKIP] Passos 1 a 4 ignorados automaticamente.")
        else:
            print(f"ℹ️  Modo FastQ Tradicional detectado para {id_amostra}.")
            
            padrao_r1 = os.path.join(PASTA_BRUTOS, f"{id_amostra}_R1_*.fastq.gz")
            busca_r1 = glob.glob(padrao_r1)
            
            if not busca_r1:
                print(f"⚠️ Alerta: FastQ R1 não encontrado para a amostra {id_amostra} na pasta '{PASTA_BRUTOS}'. Pulando.")
                continue
                
            r1 = sorted(busca_r1)[0]
            r2 = r1.replace("_R1_", "_R2_")
            
            if not os.path.exists(r2):
                print(f"⚠️ Alerta: O par R2 para o arquivo {r1} não foi localizado. Pulando amostra.")
                continue

            # --- PASSO 1: Triagem de Integridade e FastQC Inicial ---
            fastqc_html_marcador = os.path.join("fastqc_results", f"{os.path.basename(r1).replace('.fastq.gz', '_fastqc.html')}")
            if not os.path.exists(fastqc_html_marcador):
                print(f"📊 [PASSO 1] Executando triagem e controle de qualidade inicial...")
                triagem.pipeline_par(r1, r2)
            else:
                print(f"⏩ [SKIP] Passo 1 (Triagem/FastQC) já executado para {id_amostra}.")

            # --- PASSO 2: Limpeza de Adaptadores (Trimmomatic) ---
            r1_trim = os.path.join("trimmed_data", os.path.basename(r1).replace(".fastq.gz", "_trim_P.fastq.gz"))
            r2_trim = os.path.join("trimmed_data", os.path.basename(r2).replace(".fastq.gz", "_trim_P.fastq.gz"))

            if not os.path.exists(r1_trim) or not os.path.exists(r2_trim):
                print(f"✂️ [PASSO 2] Executando Trimmomatic para remoção de adaptadores...")
                limpeza.realizar_limpeza(r1, r2)
            else:
                print(f"⏩ [SKIP] Arquivos limpos (Trimmomatic) já localizados para {id_amostra}.")

            # --- PASSO 3: Alinhamento, Indexação e Mosdepth ---
            bam_gerado = os.path.join("alignment_results", f"{id_amostra}_sorted.bam")
            if not os.path.exists(bam_gerado):
                print(f"🧬 [PASSO 3] Executando Alinhamento (BWA-MEM) e Mosdepth...")
                alinhamento.realizar_alinhamento(r1_trim, r2_trim, REFERENCIA)
            else:
                print(f"⏩ [SKIP] Arquivo BAM alinhado e indexado já localizado para {id_amostra}.")

            # --- PASSO 4: Chamada de Variantes Somáticas (Mutect2) ---
            print(f"🎯 [PASSO 4] Iniciando chamada de variantes com GATK Mutect2...")
            variantes_mutect2.rodar_mutect2_smart_target(
                id_amostra=id_amostra,
                referencia=REFERENCIA,
                bam_entrada=bam_gerado,
                arquivo_bed=ARQUIVO_BED,
                pasta_saida=PASTA_VARIANTES
            )

        # --- PASSO 5: ANNOVAR + CancerVar (Comum para ambos os modos) ---
        PASTA_ANOTACAO_AMOSTRA = os.path.join(PASTA_ANOTACAO, f"anotacao_{id_amostra}")
        cancervar_final_txt = os.path.join(PASTA_ANOTACAO_AMOSTRA, f"{id_amostra}_cancervar.output.hg38_multianno.txt.cancervar")
        
        if not os.path.exists(cancervar_final_txt):
            print(f"🏷️ [PASSO 5] Iniciando anotação e predição ANNOVAR + CancerVar...")
            print(f"🔍 Usando como entrada: {os.path.basename(vcf_entrada_real)}")
            anotacao_variantes.rodar_anotacao(
                id_amostra=id_amostra,
                vcf_entrada=vcf_entrada_real,
                pasta_saida=PASTA_ANOTACAO_AMOSTRA,
                cancervar_py=CANCERVAR_PY,
                cancervar_config=CANCERVAR_CONFIG
            )
        else:
            print(f"⏩ [SKIP] Anotação CancerVar já existente em {PASTA_ANOTACAO_AMOSTRA}.")

        # --- PASSO 6: Relatórios Clínicos Finais (Excel/CSV) ---
        PASTA_RELATORIOS = "06_relatorios_finais"
        pasta_destino_amostra = os.path.join(PASTA_RELATORIOS, id_amostra)
        relatorio_excel_marcador = os.path.join(pasta_destino_amostra, f"{id_amostra}_relatorio_final.xlsx")

        if not os.path.exists(relatorio_excel_marcador):
            if os.path.exists(cancervar_final_txt):
                print(f"📊 [PASSO 6] Formatando e movendo tabelas para a pasta de relatórios...")
                pasta_destino_amostra_abs = os.path.abspath(pasta_destino_amostra)
                cancervar_final_txt_abs = os.path.abspath(cancervar_final_txt)
                script_conversor_abs = os.path.abspath(SCRIPT_CONVERSOR)
                
                comando_conversao = f"python3 '{script_conversor_abs}' '{cancervar_final_txt_abs}' '{pasta_destino_amostra_abs}' '{id_amostra}'"
                executar_script(comando_conversao)
            else:
                print(f"⚠️ Alerta: Arquivo final .cancervar não localizado para conversão em {cancervar_final_txt}.")
        else:
            print(f"⏩ [SKIP] Relatório final em Excel já gerado para {id_amostra}.")
            
        tempo_de_corrida_amostra = time.time() - tempo_inicio_amostra
        tempo_parcial_global = time.time() - tempo_inicio_global
        print(f"🏁 Amostra {id_amostra} processada nesta etapa em: {tempo_de_corrida_amostra:.2f} segundos.")
        print(f"⏱️ Tempo total acumulado do pipeline: {tempo_parcial_global:.2f} segundos.")

    tempo_total_global = time.time() - tempo_inicio_global
    print("\n======================================================================")
    print(f"🎉 PIPELINE CONCLUÍDO COM SUCESSO! Tempo Total: {tempo_total_global:.2f} segundos.")
    print("======================================================================")
