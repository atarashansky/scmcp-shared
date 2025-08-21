import os
import json
from typing import Optional, List, Dict
from fastmcp import FastMCP
from pydantic import Field
import pandas as pd
from anndata import AnnData

util_mcp = FastMCP("Util-Server")


def get_flat_dir_structure(
    path: str, ignore_hidden: bool = True, file_types: Optional[List[str]] = None
) -> Dict[str, Dict[str, List[str]]]:
    structure = {}
    for root, dirs, files in os.walk(path):
        rel_root = os.path.relpath(root, path)
        if rel_root == ".":
            rel_root = ""
        # 过滤隐藏文件和目录
        dirs[:] = [d for d in dirs if not (ignore_hidden and d.startswith("."))]
        if file_types:
            files = [f for f in files if any(f.endswith(ext) for ext in file_types)]
        files = [f for f in files if not (ignore_hidden and f.startswith("."))]
        structure[rel_root] = {"dirs": dirs, "files": files}
    return structure


def get_path_info(
    path: str, ignore_hidden: bool = True, file_types: Optional[List[str]] = None
) -> str:
    if not os.path.exists(path):
        return json.dumps({"error": "路径不存在"}, ensure_ascii=False)
    if os.path.isfile(path):
        if file_types and not any(path.endswith(ext) for ext in file_types):
            return json.dumps({"error": "文件类型不匹配"}, ensure_ascii=False)
        return json.dumps({"type": "file", "path": path}, ensure_ascii=False)
    elif os.path.isdir(path):
        structure = get_flat_dir_structure(path, ignore_hidden, file_types)
        return json.dumps(
            {"type": "directory", "path": path, "structure": structure},
            ensure_ascii=False,
            indent=2,
        )
    else:
        return json.dumps({"error": "未知类型"}, ensure_ascii=False)


@util_mcp.tool(tags={"util"})
def get_path_structure(
    path: str = Field(description="The path to get the structure of"),
) -> str:
    """get the directory structure of a path"""
    return get_path_info(path)


@util_mcp.tool(tags={"deseq2", "differential_expression"})
def run_deseq2_differential_expression_tool(
    adata: AnnData = Field(description="AnnData object with count data"),
    contrast: List[str] = Field(description="Contrast specification as [column_name, group1, group2]. group1 vs group2 comparison will be performed."),
    alpha: float = Field(default=0.05, description="Significance threshold for adjusted p-values. Must be between 0 and 1.", ge=0.0, le=1.0),
    lfc_threshold: float = Field(default=0.0, description="Log fold change threshold for significance. Must be >= 0.", ge=0.0),
    n_cpus: int = Field(default=1, description="Number of CPUs to use for analysis.", ge=1),
    layer: Optional[str] = Field(default=None, description="Layer of AnnData to use for count data. If None, uses adata.X. Can specify layer name from adata.layers.")
) -> pd.DataFrame:
    """
    Run pyDeSeq2 differential expression analysis on AnnData object.
    
    This function performs differential gene expression analysis using pyDeSeq2 
    and returns a table of statistically significant differentially expressed genes.
    
    Args:
        adata: AnnData object containing raw count data (non-negative integers)
        contrast: List of [condition_column, group1, group2] specifying the comparison
        alpha: Significance threshold for adjusted p-values (default: 0.05)
        lfc_threshold: Minimum absolute log2 fold change threshold (default: 0.0)
        n_cpus: Number of CPU cores to use for parallel processing (default: 1)
        layer: Layer of AnnData to use for count data. If None, uses adata.X (default: None)
        
    Returns:
        pandas.DataFrame: Table of significantly differentially expressed genes with columns:
            - gene_name: Gene identifier
            - log2FoldChange: Log2 fold change (group1 vs group2)
            - pvalue: Raw p-value
            - padj: Benjamini-Hochberg adjusted p-value
            - contrast: Contrast description (group1_vs_group2)
            - condition_column: Name of the condition column used
            
    Raises:
        ImportError: If pydeseq2 is not installed
        ValueError: If input validation fails (missing columns, groups, negative values)
    """
    from ..util import run_deseq2_differential_expression
    
    return run_deseq2_differential_expression(
        adata=adata,
        contrast=contrast,
        alpha=alpha,
        lfc_threshold=lfc_threshold,
        n_cpus=n_cpus,
        layer=layer
    )
