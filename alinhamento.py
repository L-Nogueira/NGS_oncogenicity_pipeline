import os
import sys
import subprocess

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
    # -t 4: usa 4 threads | -M: flag de compatibilidade para Picard/GATK
    print("--- [1/3] Mapeando reads contra o genoma (BWA-MEM) ---")
    comando_bwa = f"bwa mem -t 4 -M -R '@RG\\tID:{id_amostra}\\tSM:{id_amostra}\\tPL:ILLUMINA' {referencia} {r1} {r2} | samtools view ..."
    rodar_comando(cmd_bwa)

    # 2. Conversão para BAM e Ordenação
    # Arquivos BAM são binários (menores) e o 'sort' organiza por posição no cromossomo
    print("--- [2/3] Convertendo para BAM e ordenando (Samtools) ---")
    cmd_sort = f"samtools sort -@ 4 -o {bam_sorted} {sam_file}"
    rodar_comando(cmd_sort)

    # 3. Indexação do BAM final
    # Cria o arquivo .bai necessário para visualizar no IGV
    print("--- [3/3] Criando índice do arquivo alinhado (.bai) ---")
    cmd_index = f"samtools index {bam_sorted}"
    rodar_comando(cmd_index)

    # Limpeza do arquivo SAM (que é muito grande)
    if os.path.exists(bam_sorted):
        os.remove(sam_file)
        print(f"\n✅ Sucesso! O alinhamento ordenado está em: {bam_sorted}")

if __name__ == "__main__":
    if len(sys.argv) == 4:
        realizar_alinhamento(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("Uso: python3 alinhamento.py <R1_trim.fastq.gz> <R2_trim.fastq.gz> <referencia.fa>")
