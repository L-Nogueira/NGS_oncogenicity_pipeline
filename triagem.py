import gzip
import hashlib
import sys
import os

def calcular_md5(caminho):
    hasher = hashlib.md5()
    with open(caminho, 'rb') as f:
        for bloco in iter(lambda: f.read(4096), b""):
            hasher.update(bloco)
    return hasher.hexdigest()

def contar_reads(caminho):
    with gzip.open(caminho, 'rt') as f:
        for i, linha in enumerate(f):
            pass
        return (i + 1) // 4

def pipeline_par(r1, r2):
    print(f"\n=== Iniciando Processamento de Par Real ===")

    # Processando R1
    print(f"Lendo R1: {r1}...")
    md5_r1 = calcular_md5(r1)
    reads_r1 = contar_reads(r1)

    # Processando R2
    print(f"Lendo R2: {r2}...")
    md5_r2 = calcular_md5(r2)
    reads_r2 = contar_reads(r2)

    print("-" * 30)
    print(f"RESULTADOS:")
    print(f"R1 MD5: {md5_r1} | Reads: {reads_r1}")
    print(f"R2 MD5: {md5_r2} | Reads: {reads_r2}")

    # Verificação de Integridade de Par
    if reads_r1 == reads_r2:
        print("\n✅ SUCESSO: O número de reads em R1 e R2 é idêntico.")

        print("\n[3/3] Iniciando FastQC (Controle de Qualidade)...")

        pasta_output = "fastqc_results"

        if not os.path.exists(pasta_output):
                os.makedirs(pasta_output)
                print(f"Diretório '{pasta_output}' criado.")

        os.system(f"fastqc {r1} {r2} -o {pasta_output}")
        print(f"\nFastQC finalizado! Verifique os arquivos .html na pasta {pasta_output}/")
    else:
        print("\n❌ ERRO: R1 e R2 possuem números diferentes de reads!")

if __name__ == "__main__":
    if len(sys.argv) == 3:
        pipeline_par(sys.argv[1], sys.argv[2])
    else:
        print("Uso: python3 triagem.py amostra_R1.fastq.gz amostra_R2.fastq.gz")
