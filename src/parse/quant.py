from src.parse.flashquant import parseFLASHQuantOutput
from src.parse.flashquant import connectTraceWithResult


def parseQuant(quant_tsv, trace_tsv, conflict_tsv=None):
        quant_df, trace_df, resolution_df = parseFLASHQuantOutput(
            quant_tsv, trace_tsv, conflict_tsv
        )
        quant_df = connectTraceWithResult(quant_df, trace_df)
        
        results = {'quant_dfs' : quant_df}
        if resolution_df is not None:
            results['conflict_resolution_dfs'] = resolution_df
        
        return results
