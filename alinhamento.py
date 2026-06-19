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

    # 2. Conversão para BAM
    print("--- [2/3] Convertendo SAM para BAM ---")
    bam_file = sam_file.replace(".sam", ".bam")
    cmd_view = f"samtools view -bS '{sam_file}' > '{bam_file}'"
    rodar_comando(cmd_view)

    # 3. Ordenação do BAM e Indexação
    print("--- [3/3] Ordenando e indexando arquivo BAM (Samtools) ---")
    cmd_sort = f"samtools sort -@ 4 '{bam_file}' -o '{bam_sorted}'"
    rodar_comando(cmd_sort)
    
    cmd_index = f"samtools index '{bam_sorted}'"
    rodar_comando(cmd_index)

    # Remoção de arquivos intermediários volumosos para poupar espaço
    if os.path.exists(sam_file): os.remove(sam_file)
    if os.path.exists(bam_file): os.remove(bam_file)
    print(f"✅ Alinhamento concluído! BAM indexado criado em: {bam_sorted}")

    # ==================================================================
    # REVISÃO SMART TARGET: Identificação Automática do arquivo BED
    # ==================================================================
    pasta_bed = "/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/arquivo .bed/"
    arquivos_bed = glob.glob(os.path.join(pasta_bed, "*.bed"))
    
    # Filtrar para ignorar o arquivo temporário FIXED caso ele já exista na pasta
    arquivos_bed = [b for b in arquivos_bed if "targeted_regions_FIXED.bed" not in os.path.basename(b)]

    if not arquivos_bed:
        print(f"❌ Erro crítico: Nenhum arquivo .bed localizado na pasta {pasta_bed}")
        sys.exit(1)
        
    # Seleção automática: adota o primeiro arquivo BED encontrado na pasta
    bed_selecionado = arquivos_bed[0]
    print(f"🎯 Arquivo BED selecionado automaticamente: {os.path.basename(bed_selecionado)}")
    
    bed_corrigido = os.path.join(pasta_bed, "targeted_regions_FIXED.bed")
    print("🔧 Padronizando nomenclatura do arquivo BED (removendo prefixo 'chr')...")
    
    cmd_sed = f"sed 's/^chr//' '{bed_selecionado}' > '{bed_corrigido}'"
    subprocess.run(cmd_sed, shell=True, check=True)

    # ------------------------------------------------------------------
    # REDIRECIONAMENTO DO MOSDEPTH PARA A PASTA 06_RELATORIOS_FINAIS
    # ------------------------------------------------------------------
    pasta_mosdepth_out = os.path.join("06_relatorios_finais", id_amostra)
    os.makedirs(pasta_mosdepth_out, exist_ok=True)
    
    prefixo_mosdepth = os.path.abspath(os.path.join(pasta_mosdepth_out, f"{id_amostra}_focado_alvo"))
    bam_sorted_abs = os.path.abspath(bam_sorted)
    bed_corrigido_abs = os.path.abspath(bed_corrigido)
    
    print(f"🧬 Rodando Mosdepth (Salvando em {pasta_mosdepth_out})...")
    cmd_mosdepth = f"mosdepth -n --fast-mode -b '{bed_corrigido_abs}' '{prefixo_mosdepth}' '{bam_sorted_abs}'"
    
    try:
        subprocess.run(cmd_mosdepth, shell=True, check=True)
        print(f"✅ Relatório do Mosdepth gerado com sucesso!")
        print(f"📁 Arquivo de resumo disponível em: {prefixo_mosdepth}.mosdepth.summary.txt")
    except Exception as e:
        print(f"⚠️ Falha não impeditiva no Mosdepth: {e}")

    return bam_sorted
