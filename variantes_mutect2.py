#!/usr/bin/env python3
import os
import subprocess
import sys

def rodar_mutect2_smart_target(id_amostra, referencia, bam_entrada, arquivo_bed, pasta_saida):
    """
    Roda a chamada de variantes somáticas Tumor-Only usando o GATK Mutect2
    restringindo o processamento às regiões do arquivo BED padronizado (Smart Target).
    """
    os.makedirs(pasta_saida, exist_ok=True)
    
    vcf_bruto = os.path.join(pasta_saida, f"{id_amostra}_unfiltered.vcf")
    vcf_filtrado = os.path.join(pasta_saida, f"{id_amostra}_variants.vcf")
    stats_mutect = f"{vcf_bruto}.stats" 
    
    print(f"🧬 [Mutect2] Iniciando chamada de variantes somáticas para: {id_amostra}")
    
    # ------------------------------------------------------------------
    # CORREÇÃO: Buscar o arquivo BED temporário corrigido (sem 'chr') gerado no alinhamento
    # Isso evita incompatibilidade de cromossomos e lida com o espaço na pasta
    # ------------------------------------------------------------------
    diretorio_bed = os.path.dirname(arquivo_bed)
    bed_corrigido = os.path.join(diretorio_bed, "targeted_regions_FIXED.bed")
    
    # Se o alinhamento já limpou ou não gerou o FIXED, usamos um fallback seguro
    if os.path.exists(bed_corrigido):
        bed_alvo = bed_corrigido
        print(f"🎯 Usando arquivo BED padronizado para o GATK: {os.path.basename(bed_alvo)}")
    else:
        bed_alvo = arquivo_bed
        print(f"⚠️ Aviso: Arquivo BED modificado não encontrado. Usando o original.")

    # Comando protegido por aspas duplas \" para caminhos com espaços no sistema de arquivos
    cmd_mutect = (
        f"gatk Mutect2 "
        f"-R \"{referencia}\" "
        f"-I \"{bam_entrada}\" "
        f"-L \"{bed_alvo}\" "
        f"-O \"{vcf_bruto}\""
    )
    
    print(f"🏃‍♂️ Executando Passo 1/2 (Mutect2)...")
    res_mutect = subprocess.run(cmd_mutect, shell=True, capture_output=True, text=True)
    
    if res_mutect.returncode != 0:
        print(f"❌ Erro no MuTect2 para a amostra {id_amostra}:", file=sys.stderr)
        print(res_mutect.stderr, file=sys.stderr)
        return False

    cmd_filter = (
        f"gatk FilterMutectCalls "
        f"-R \"{referencia}\" "
        f"-V \"{vcf_bruto}\" "
        f"-O \"{vcf_filtrado}\""
    )
    
    print(f"🧹 Executando Passo 2/2 (FilterMutectCalls)...")
    res_filter = subprocess.run(cmd_filter, shell=True, capture_output=True, text=True)
    
    if res_filter.returncode != 0:
        print(f"❌ Erro no FilterMutectCalls para a amostra {id_amostra}:", file=sys.stderr)
        print(res_filter.stderr, file=sys.stderr)
        return False
        
    print(f"♻️  Limpando arquivos intermediários da amostra {id_amostra}...")
    for arquivo_temp in [vcf_bruto, f"{vcf_bruto}.idx", stats_mutect]:
        if os.path.exists(arquivo_temp):
            os.remove(arquivo_temp)
            
    return True
