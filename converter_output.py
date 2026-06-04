import pandas as pd
import sys
import os

def converter_cancervar(input_txt, output_dir, id_amostra):
    print(f"📊 Lendo e corrigindo arquivo: {os.path.basename(input_txt)}")
    
    try:
        # Lemos o arquivo bruto separando por tabulação
        df = pd.read_csv(input_txt, sep='\t', low_memory=False)
        
        # CORREÇÃO 1: Remover a linha fantasma de cabeçalho repetido se existir
        if df.iloc[0, 1] == 'Start':
            df = df.drop(df.index[0]).reset_index(drop=True)
            
        # CORREÇÃO 2: Ajustar nomes das colunas (remover espaços)
        df.columns = [c.strip() for c in df.columns]
        
        # Ajusta a primeira coluna para ficar legível se vier com o nome longo do CancerVar
        if df.columns[0].startswith('CancerVar:'):
            df.rename(columns={df.columns[0]: 'Interpretation_Details'}, inplace=True)

    except Exception as e:
        print(f"❌ Erro ao ler e processar o arquivo: {e}")
        return

    colunas = list(df.columns)
    
    # Ordem prioritária de colunas para a bancada
    colunas_prioritarias = [
        'Interpretation_Details',
        'Chr', 'Start', 'End', 'Ref', 'Alt',
        'Gene.refGene', 'Func.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
        'Gene.ensGene', 'AAChange.ensGene',
        'clinvar: Clinvar', 'cosmic91'
    ]
    
    colunas_para_o_inicio = [c for c in colunas_prioritarias if c in colunas]
    outras_colunas = [c for c in colunas if c not in colunas_para_o_inicio]
    df_reorganizado = df[colunas_para_o_inicio + outras_colunas]

    # GARANTIR DIRETÓRIO: Cria a pasta final da amostra se ela não existir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"📁 Nova pasta criada para os relatórios: {output_dir}")

    # Define o caminho completo dos arquivos finais usando o ID da amostra
    csv_out = os.path.join(output_dir, f"{id_amostra}_relatorio_final.csv")
    xlsx_out = os.path.join(output_dir, f"{id_amostra}_relatorio_final.xlsx")

    # Salva em CSV
    df_reorganizado.to_csv(csv_out, index=False, sep=';')
    print(f"✅ CSV estruturado gerado: {csv_out}")

    # Salva em XLSX (Excel)
    with pd.ExcelWriter(xlsx_out, engine='openpyxl') as writer:
        df_reorganizado.to_excel(writer, index=False, sheet_name='Variantes_Somáticas')
    print(f"✅ Planilha Excel (.xlsx) gerada: {xlsx_out}")

if __name__ == "__main__":
    # Mantém valores padrão compatíveis caso seja executado manualmente
    input_padrao = "/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/05_anotacao/1181_S12_L001_cancervar.output.hg38_multianno.txt.cancervar"
    dir_out_padrao = "/home/l.nogueira/laboratorio_bioinfo/projetos_miseq_real/06_relatorios_finais/1181_S12_L001"
    id_padrao = "1181_S12_L001"
    
    path_in = sys.argv[1] if len(sys.argv) > 1 else input_padrao
    path_out_dir = sys.argv[2] if len(sys.argv) > 2 else dir_out_padrao
    sample_id = sys.argv[3] if len(sys.argv) > 3 else id_padrao

    if os.path.exists(path_in):
        converter_cancervar(path_in, path_out_dir, sample_id)
    else:
        print(f"❌ Arquivo de entrada não encontrado: {path_in}")
