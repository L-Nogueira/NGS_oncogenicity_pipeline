import os
import sys
import subprocess
import glob

def rodar_comando(comando):
    """Função auxiliar para executar comandos no terminal e tratar erros."""
    try:
        subprocess.run(comando, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Erro ao executar: {comando}\nErro: {e}")
        sys.exit(1)

def realizar_alinhamento(r1, r2, referencia):
    print(f"\n🚀 Iniciando Alinhamento: {os.path.basename(r1)}")
    
    # Criar pasta de resultados se não existir
    pasta_out = "alignment_results"
    if not os.path.exists(pasta_out):
        os.makedirs(pasta_out)

    # Definir nomes dos arquivos baseados na amostra
    id_amostra = os.path.basename(r1).split('_R1')[0]
    sam_file = os.path.join(pasta_out, f"{id_amostra}.sam")
    bam_sorted = os.path.join(pasta_out, f"{id_amostra}_sorted.bam")

    # 1. Alinhamento com BWA-MEM
    print("--- [1/3] Mapeando reads contra o genoma (BWA-MEM) ---")
    cmd_bwa = f"bwa mem -t 4 -M -R '@RG\\tID:{id_amostra}\\tSM:{id_amostra}\\tPL:ILLUMINA' '{referencia}' '{r1}' '{r2}' > '{sam_file}'"
    rodar_comando(cmd_bwa)

    # 2. Conversão SAM para BAM, Ordenação e Indexação (Samtools)
    print("--- [2/3] Convertendo, ordenando e indexando (Samtools) ---")
    cmd_samtools = f"samtools sort -@ 4 -o '{bam_sorted}' '{sam_file}' && samtools index '{bam_sorted}'"
    rodar_comando(cmd_samtools)

    # Remover o arquivo SAM residual para economizar espaço em disco
    if os.path.exists(sam_file):
        os.remove(sam_file)

    # 3. Localização automática e correção do arquivo BED
    print("--- [3/3] Processando Controle de Qualidade de Cobertura ---")
    pasta_bed = "/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/arquivo .bed"
    arquivos_bed = glob.glob(os.path.join(pasta_bed, "*.bed"))
    
    if not arquivos_bed:
        print("❌ Erro: Nenhum arquivo .bed encontrado na pasta especificada.")
        sys.exit(1)
        
    bed_selecionado = arquivos_bed[0]
    print(f"🎯 Arquivo BED selecionado automaticamente: {os.path.basename(bed_selecionado)}")
    
    bed_corrigido = os.path.join(pasta_bed, "targeted_regions_FIXED.bed")
    print("🔧 Padronizando nomenclatura do arquivo BED (removendo prefixo 'chr')...")
    
    cmd_sed = f"sed 's/^chr//' '{bed_selecionado}' > '{bed_corrigido}'"
    subprocess.run(cmd_sed, shell=True, check=True)

    # ------------------------------------------------------------------
    # CHAMADA PROTEGIDA DO MOSDEPTH VIA TRY-EXCEPT
    # ------------------------------------------------------------------
    pasta_mosdepth_out = os.path.join("06_relatorios_finais", id_amostra)
    os.makedirs(pasta_mosdepth_out, exist_ok=True)
    
    prefixo_mosdepth = os.path.abspath(os.path.join(pasta_mosdepth_out, f"{id_amostra}_focado_alvo"))
    bam_sorted_abs = os.path.abspath(bam_sorted)
    bed_corrigido_abs = os.path.abspath(bed_corrigido)
    
    if os.path.exists(bam_sorted_abs):
        try:
            print(f"🧬 Rodando Mosdepth (Salvando em {pasta_mosdepth_out})...")
            cmd_mosdepth = f"mosdepth -n --fast-mode -b '{bed_corrigido_abs}' '{prefixo_mosdepth}' '{bam_sorted_abs}'"
            
            res_mosdepth = subprocess.run(cmd_mosdepth, shell=True, capture_output=True, text=True)
            if res_mosdepth.returncode == 0:
                print("📊 [SUCESSO] Relatório do Mosdepth gerado perfeitamente!")
            else:
                print(f"⚠️ Erro interno do Mosdepth (Execução continuará):\n{res_mosdepth.stderr}")
        except Exception as e:
            print(f"⚠️ Falha crítica ao disparar o Mosdepth (Execução continuará): {e}")
    else:
        print(f"❌ Erro: Arquivo BAM indexado não encontrado em {bam_sorted_abs}. Mosdepth ignorado.")

    print(f"🏁 Alinhamento e Controle de Qualidade concluídos para a amostra {id_amostra}.\n")
