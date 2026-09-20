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

def realizar_variant_calling(bam_original, bed_original, referencia):
    print(f"\n🧬 Iniciando Variant Calling para: {os.path.basename(bam_original)}")
    
    # 1. Definição e criação das pastas de saída
    pasta_03 = "03_alinhamento"
    pasta_04 = "04_variantes"
    for pasta in [pasta_03, pasta_04]:
        if not os.path.exists(pasta):
            os.makedirs(pasta)

    # Extrai o ID da amostra (ex: 1170) baseado no nome do arquivo BAM
    id_amostra = os.path.basename(bam_original).split('_')[0]
    
    # Definição dos caminhos dos arquivos gerados
    bed_formatado = os.path.join(pasta_03, "targeted_regions_ensembl.bed")
    bam_targets = os.path.join(pasta_03, f"{id_amostra}_targets.bam")
    bcf_saida = os.path.join(pasta_04, f"{id_amostra}_variants.bcf")
    vcf_final = os.path.join(pasta_04, f"{id_amostra}_variants.vcf")

    # 2. Normalização do BED (Removendo 'chr' para bater com o padrão do BAM)
    print("--- [1/5] Normalizando nomenclatura do BED (Estilo Ensembl) ---")
    cmd_sed = f"sed 's/^chr//' '{bed_original}' > '{bed_formatado}'"
    rodar_comando(cmd_sed)

    # 3. Filtrar o BAM original usando as regiões do BED
    print("--- [2/5] Filtrando arquivo BAM com base nas regiões do BED ---")
    cmd_view = f"samtools view -b -L '{bed_formatado}' '{bam_original}' > '{bam_targets}'"
    rodar_comando(cmd_view)

    # 4. Indexar o BAM filtrado
    print("--- [3/5] Indexando o arquivo BAM filtrado (.bai) ---")
    cmd_index = f"samtools index '{bam_targets}'"
    rodar_comando(cmd_index)

    # 5. Chamada de Variantes (mpileup + call)
    print("--- [4/5] Executando bcftools mpileup e call ---")
    cmd_bcftools = (
        f"bcftools mpileup -f '{referencia}' -R '{bed_formatado}' '{bam_targets}' 2>/dev/null | "
        f"bcftools call -mv -Ob -o '{bcf_saida}'"
    )
    rodar_comando(cmd_bcftools)

    # 6. Conversão de BCF para VCF legível
    print("--- [5/5] Convertendo resultado binário para VCF plano ---")
    cmd_view_vcf = f"bcftools view '{bcf_saida}' > '{vcf_final}'"
    rodar_comando(cmd_view_vcf)

    # 7. Limpeza de arquivos temporários e intermediários (Boas práticas de ADS)
    print("--- Limpando arquivos temporários para economizar espaço ---")
    if os.path.exists(bed_formatado): os.remove(bed_formatado)
    if os.path.exists(bcf_saida): os.remove(bcf_saida)
    if os.path.exists(bam_targets): os.remove(bam_targets)
    if os.path.exists(f"{bam_targets}.bai"): os.remove(f"{bam_targets}.bai")

    print(f"✅ SUCESSO: O arquivo final foi gerado em: {vcf_final}\n")

if __name__ == "__main__":
    # Verifica se os argumentos necessários foram passados
    if len(sys.argv) < 4:
        print("Uso correto: python3 variantes.py <arquivo_BAM> <arquivo_BED> <referencia_FASTA>")
        sys.exit(1)
        
    bam_in = sys.argv[1]
    bed_in = sys.argv[2]
    ref_in = sys.argv[3]
    
    realizar_variant_calling(bam_in, bed_in, ref_in)
