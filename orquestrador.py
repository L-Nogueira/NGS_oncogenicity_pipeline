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
# --- CONFIGURAÇÃO DE CAMINHOS E EXECUTÁVEIS ---
# ======================================================================
PASTA_BRUTOS = "dados_brutos"
REFERENCIA = os.path.expanduser("~/laboratorio_bioinfo/genomas_referencia/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
PASTA_VARIANTES = "04_variantes"
PASTA_ANOTACAO = "05_anotacao" 

# Arquivos e caminhos de suporte - Alterado para compatibilidade nativa com a automação
ARQUIVO_BED = "/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/arquivo .bed/targeted_regions_FIXED.bed"
ANNOVAR_DIR = "/home/l.nogueira/laboratorio_bioinfo/softwares/annovar"

# CANCERVAR ATIVO
CANCERVAR_DIR = "/home/l.nogueira/laboratorio_bioinfo/softwares/CancerVar"
CANCERVAR_PY = os.path.join(CANCERVAR_DIR, "CancerVar.py")
CANCERVAR_CONFIG = os.path.join(CANCERVAR_DIR, "config.ini")

# Caminho do script de conversão amigável (Excel/CSV)
SCRIPT_CONVERSOR = "/home/l.nogueira/laboratorio_bioinfo/scripts_bioinfo/converter_output.py"

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

    # Coleta de arquivos R1 de forma robusta
    arquivos_r1 = sorted(glob.glob(os.path.join(PASTA_BRUTOS, "*_R1_*.fastq.gz")))

    if not arquivos_r1:
        print(f"❌ Erro: Nenhum arquivo pareado R1 encontrado na pasta '{PASTA_BRUTOS}'.")
        sys.exit(1)

    print(f"📦 Total de amostras identificadas para processamento: {len(arquivos_r1)}")

    for r1 in arquivos_r1:
        tempo_inicio_amostra = time.time()
        
        # Identificar par R2 correspondente mudando estritamente a ocorrência de _R1_ para _R2_
        r2 = r1.replace("_R1_", "_R2_")
        if not os.path.exists(r2):
            print(f"⚠️ Alerta: O par R2 para o arquivo {r1} não foi localizado. Pulando amostra.")
            continue

        id_amostra = obter_id_amostra(r1)
        print(f"\n{'-'*80}")
        print(f"🔬 PROCESSANDO AMOSTRA: {id_amostra}")
        print(f"{'-'*80}")

        # --- PASSO 1: Triagem de Integridade e FastQC Inicial (Com SKIP inteligente) ---
        # Define um arquivo alvo gerado pelo FastQC para servir de marcador de sucesso
        fastqc_html_marcador = os.path.join("fastqc_results", f"{os.path.basename(r1).replace('.fastq.gz', '_fastqc.html')}")
        
        if not os.path.exists(fastqc_html_marcador):
            print(f"📊 [PASSO 1] Executando triagem e controle de qualidade inicial...")
            triagem.pipeline_par(r1, r2)
        else:
            print(f"⏩ [SKIP] Passo 1 (Triagem/FastQC) já executado para {id_amostra}.")

        # --- PASSO 2: Limpeza de Adaptadores (Trimmomatic) (Com SKIP Inteligente) ---
        r1_trim = os.path.join("trimmed_data", os.path.basename(r1).replace(".fastq.gz", "_trim_P.fastq.gz"))
        r2_trim = os.path.join("trimmed_data", os.path.basename(r2).replace(".fastq.gz", "_trim_P.fastq.gz"))

        if not os.path.exists(r1_trim) or not os.path.exists(r2_trim):
            print(f"✂️ [PASSO 2] Executando Trimmomatic para remoção de adaptadores...")
            limpeza.realizar_limpeza(r1, r2)
        else:
            print(f"⏩ [SKIP] Arquivos limpos (Trimmomatic) já localizados para {id_amostra}.")

        # --- PASSO 3: Alinhamento, Indexação e Mosdepth (Com SKIP Inteligente) ---
        # IMPORTANTE: Daqui para frente usamos r1_trim e r2_trim!
        bam_gerado = os.path.join("alignment_results", f"{id_amostra}_sorted.bam")
        
        if not os.path.exists(bam_gerado):
            print(f"🧬 [PASSO 3] Executando Alinhamento (BWA-MEM) e Mosdepth...")
            # Corrigido internamente em alinhamento.py com \\t
            alinhamento.realizar_alinhamento(r1_trim, r2_trim, REFERENCIA)
        else:
            print(f"⏩ [SKIP] Arquivo BAM alinhado e indexado já localizado para {id_amostra}.")

        # --- PASSO 4: Chamada de Variantes Somáticas (Mutect2) (Com SKIP Inteligente) ---
        vcf_gerado = os.path.join(PASTA_VARIANTES, f"{id_amostra}_variants.vcf")
        
        if not os.path.exists(vcf_gerado):
            print(f"🎯 [PASSO 4] Iniciando chamada de variantes com GATK Mutect2...")
            variantes_mutect2.rodar_mutect2_smart_target(
                id_amostra=id_amostra,
                referencia=REFERENCIA,
                bam_entrada=bam_gerado,
                arquivo_bed=ARQUIVO_BED,
                pasta_saida=PASTA_VARIANTES
            )
        else:
            print(f"⏩ [SKIP] Chamada de variantes (VCF) já localizada para {id_amostra}.")

        # --- PASSO 5: ANNOVAR + CancerVar (Com SKIP Inteligente) ---
        cancervar_final_txt = os.path.join(PASTA_ANOTACAO, f"{id_amostra}_cancervar.output.hg38_multianno.txt.cancervar")
        
        if not os.path.exists(cancervar_final_txt):
            print(f"🏷️ [PASSO 5] Iniciando anotação e predição ANNOVAR + CancerVar...")
            anotacao_variantes.rodar_anotacao(
                id_amostra=id_amostra,
                vcf_entrada=vcf_gerado,
                pasta_saida=PASTA_ANOTACAO,
                cancervar_py=CANCERVAR_PY,
                cancervar_config=CANCERVAR_CONFIG
            )
        else:
            print(f"⏩ [SKIP] Anotação CancerVar já existente para {id_amostra}.")

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
                
                comando_conversao = f"python3 \"{script_conversor_abs}\" \"{cancervar_final_txt_abs}\" \"{pasta_destino_amostra_abs}\" \"{id_amostra}\""
                executar_script(comando_conversao)
            else:
                print(f"⚠️ Alerta: Arquivo final .cancervar não localizado para conversão. Verifique logs internos.")
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
