import os
import sys

def realizar_limpeza(r1, r2):
    print(f"\n=== Iniciando Etapa de Limpeza (Trimming) ===")
    
    # Definição das pastas de saída
    pasta_data = "trimmed_data"
    pasta_qc_pos = "fastqc_pos"

    # Criando as pastas se não existirem
    for pasta in [pasta_data, pasta_qc_pos]:
        if not os.path.exists(pasta):
            os.makedirs(pasta)
            print(f"Diretório '{pasta}' criado.")

    # Definindo caminhos dos arquivos de saída dentro da pasta 'trimmed_data'
    # Usamos o os.path.basename para não carregar o caminho antigo no nome do arquivo
    nome_r1 = os.path.basename(r1)
    nome_r2 = os.path.basename(r2)

    r1_trim_p = os.path.join(pasta_data, nome_r1.replace(".fastq.gz", "_trim_P.fastq.gz"))
    r1_trim_u = os.path.join(pasta_data, nome_r1.replace(".fastq.gz", "_trim_U.fastq.gz"))
    r2_trim_p = os.path.join(pasta_data, nome_r2.replace(".fastq.gz", "_trim_P.fastq.gz"))
    r2_trim_u = os.path.join(pasta_data, nome_r2.replace(".fastq.gz", "_trim_U.fastq.gz"))

    # Comando Trimmomatic
    # O arquivo TruSeq3-PE.fa geralmente está no path do conda, 
    # mas se der erro, precisaremos do caminho completo.
    comando_trim = (
        f"trimmomatic PE -phred33 {r1} {r2} "
        f"{r1_trim_p} {r1_trim_u} {r2_trim_p} {r2_trim_u} "
        f"ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 "
        f"LEADING:3 TRAILING:3 "
        f"SLIDINGWINDOW:4:15 "
        f"MINLEN:36"
    )

    print(f"\n[1/2] Executando Trimmomatic...")
    os.system(comando_trim)

    # Gerando o FastQC dos dados LIMPOS na pasta 'fastqc_pos'
    print(f"\n[2/2] Gerando FastQC pós-limpeza em '{pasta_qc_pos}'...")
    os.system(f"fastqc {r1_trim_p} {r2_trim_p} -o {pasta_qc_pos}")
    
    print("\n" + "="*30)
    print(f"✅ PROCESSO CONCLUÍDO")
    print(f"Output de sequências: ./{pasta_data}/")
    print(f"Relatórios de qualidade: ./{pasta_qc_pos}/")
    print("="*30)

if __name__ == "__main__":
    if len(sys.argv) == 3:
        realizar_limpeza(sys.argv[1], sys.argv[2])
    else:
        print("Uso: python3 limpeza.py R1.fastq.gz R2.fastq.gz")
        
