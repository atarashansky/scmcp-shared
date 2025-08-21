from typing import Union
import pandas as pd
from fastmcp.exceptions import ToolError
from fastmcp.tools.tool import Tool
from scmcp_shared.schema.preset.util import *
from scmcp_shared.schema.util import DeSeq2DEParam
from scmcp_shared.schema.preset import AdataInfo
from scmcp_shared.util import (
    forward_request,
    get_ads,
    add_op_log,
    deserialize_mcp_param,
)
from scmcp_shared.mcp_base import BaseMCP


class ScanpyUtilMCP(BaseMCP):
    def __init__(
        self,
        include_tools: list = None,
        exclude_tools: list = None,
        AdataInfo=AdataInfo,
    ):
        """
        Initialize ScanpyUtilMCP with optional tool filtering.

        Args:
            include_tools (list, optional): List of tool names to include. If None, all tools are included.
            exclude_tools (list, optional): List of tool names to exclude. If None, no tools are excluded.
            AdataInfo: The AdataInfo class to use for type annotations.
        """
        super().__init__("SCMCP-Util-Server", include_tools, exclude_tools, AdataInfo)

    def _tool_query_op_log(self):
        def _query_op_log(
            request: Union[QueryOpLogParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Query the adata operation log"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, QueryOpLogParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            adata = get_ads().get_adata(adinfo=adinfo)
            op_dic = adata.uns["operation"]["op"]
            opids = adata.uns["operation"]["opid"][-request.n :]
            op_list = []
            for opid in opids:
                op_list.append(op_dic[opid])
            return op_list

        return Tool.from_function(
            _query_op_log, name="query_op_log", enabled=True, tags=["preset"]
        )

    def _tool_mark_var(self):
        def _mark_var(
            request: Union[MarkVarParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """
            Determine if each gene meets specific conditions and store results in adata.var as boolean values.
            For example: mitochondrion genes startswith MT-.
            The tool should be called first when calculate quality control metrics for mitochondrion, ribosomal, harhemoglobin genes, or other qc_vars.
            """
            # Deserialize parameters
            request = deserialize_mcp_param(request, MarkVarParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("ul_mark_var", request, adinfo)
                if result is not None:
                    return result
                adata = get_ads().get_adata(adinfo=adinfo)
                var_name = request.var_name
                gene_class = request.gene_class
                pattern_type = request.pattern_type
                patterns = request.patterns
                if gene_class is not None:
                    if gene_class == "mitochondrion":
                        adata.var["mt"] = adata.var_names.str.startswith(
                            ("MT-", "Mt", "mt-")
                        )
                        var_name = "mt"
                    elif gene_class == "ribosomal":
                        adata.var["ribo"] = adata.var_names.str.startswith(
                            ("RPS", "RPL", "Rps", "Rpl")
                        )
                        var_name = "ribo"
                    elif gene_class == "hemoglobin":
                        adata.var["hb"] = adata.var_names.str.contains(
                            "^HB[^(P)]", case=False
                        )
                        var_name = "hb"
                elif pattern_type is not None and patterns is not None:
                    if pattern_type == "startswith":
                        adata.var[var_name] = adata.var_names.str.startswith(patterns)
                    elif pattern_type == "endswith":
                        adata.var[var_name] = adata.var_names.str.endswith(patterns)
                    elif pattern_type == "contains":
                        adata.var[var_name] = adata.var_names.str.contains(patterns)
                    else:
                        raise ValueError(
                            f"Did not support pattern_type: {pattern_type}"
                        )
                else:
                    raise ValueError("Please provide validated parameter")

                res = {
                    var_name: adata.var[var_name].value_counts().to_dict(),
                    "msg": f"add '{var_name}' column in adata.var",
                }
                func_kwargs = {
                    "var_name": var_name,
                    "gene_class": gene_class,
                    "pattern_type": pattern_type,
                    "patterns": patterns,
                }
                add_op_log(adata, "mark_var", func_kwargs, adinfo)
                return res
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _mark_var, name="mark_var", enabled=True, tags=["preset"]
        )

    def _tool_list_var(self):
        def _list_var(
            request: Union[ListVarParam, str, dict] = None,
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """List key columns in adata.var. It should be called for checking when other tools need var key column names as input."""
            # Deserialize parameters
            request = deserialize_mcp_param(request, ListVarParam, ListVarParam())
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("ul_list_var", request, adinfo)
                if result is not None:
                    return result
                adata = get_ads().get_adata(adinfo=adinfo)
                columns = list(adata.var.columns)
                add_op_log(adata, "list_var", {}, adinfo)
                return columns
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _list_var, name="list_var", enabled=True, tags=["preset"]
        )

    def _tool_list_obs(self):
        def _list_obs(
            request: Union[ListObsParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """List key columns in adata.obs. It should be called before other tools need obs key column names input."""
            # Deserialize parameters
            print("TEST REQUEST LOG")
            request = deserialize_mcp_param(request, ListObsParam)
            print(request)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("ul_list_obs", request, adinfo)
                if result is not None:
                    return result
                adata = get_ads().get_adata(adinfo=adinfo)
                columns = list(adata.obs.columns)
                add_op_log(adata, "list_obs", {}, adinfo)
                return columns
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _list_obs, name="list_obs", enabled=True, tags=["preset"]
        )

    def _tool_check_var(self):
        def _check_var(
            request: Union[VarNamesParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Check if genes/variables exist in adata.var_names. This tool should be called before gene expression visualizations or color by genes."""
            # Deserialize parameters
            request = deserialize_mcp_param(request, VarNamesParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("ul_check_var", request, adinfo)
                if result is not None:
                    return result
                adata = get_ads().get_adata(adinfo=adinfo)
                var_names = request.var_names
                if adata.raw is not None:
                    all_var_names = adata.raw.to_adata().var_names
                else:
                    all_var_names = adata.var_names
                result = {v: v in all_var_names for v in var_names}
                add_op_log(adata, "check_var", {"var_names": var_names}, adinfo)
                return result
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _check_var, name="check_var", enabled=True, tags=["preset"]
        )

    def _tool_merge_adata(self):
        def _merge_adata(
            request: Union[ConcatBaseParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Merge multiple adata objects."""
            # Deserialize parameters
            request = deserialize_mcp_param(request, ConcatBaseParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("ul_merge_adata", request, adinfo)
                if result is not None:
                    return result
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                kwargs = {
                    k: v for k, v in request.model_dump().items() if v is not None
                }
                merged_adata = adata.concat(ads.adata_dic, **kwargs)
                ads.adata_dic[dtype] = {}
                ads.active_id = "merged_adata"
                add_op_log(merged_adata, ad.concat, kwargs, adinfo)
                ads.adata_dic[ads.active_id] = merged_adata
                return {
                    "status": "success",
                    "message": "Successfully merged all AnnData objects",
                }
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _merge_adata, name="merge_adata", enabled=True, tags=["preset"]
        )

    def _tool_set_dpt_iroot(self):
        def _set_dpt_iroot(
            request: Union[DPTIROOTParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Set the iroot cell"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, DPTIROOTParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("ul_set_dpt_iroot", request, adinfo)
                if result is not None:
                    return result
                adata = get_ads().get_adata(adinfo=adinfo)
                diffmap_key = request.diffmap_key
                dimension = request.dimension
                direction = request.direction
                if diffmap_key not in adata.obsm:
                    raise ValueError(
                        f"Diffusion map key '{diffmap_key}' not found in adata.obsm"
                    )
                if direction == "min":
                    adata.uns["iroot"] = adata.obsm[diffmap_key][:, dimension].argmin()
                else:
                    adata.uns["iroot"] = adata.obsm[diffmap_key][:, dimension].argmax()

                func_kwargs = {
                    "diffmap_key": diffmap_key,
                    "dimension": dimension,
                    "direction": direction,
                }
                add_op_log(adata, "set_dpt_iroot", func_kwargs, adinfo)

                return {
                    "status": "success",
                    "message": f"Successfully set root cell for DPT using {direction} of dimension {dimension}",
                }
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _set_dpt_iroot, name="set_dpt_iroot", enabled=True, tags=["preset"]
        )

    def _tool_add_layer(self):
        def _add_layer(
            request: Union[AddLayerParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Add a layer to the AnnData object."""
            # Deserialize parameters
            request = deserialize_mcp_param(request, AddLayerParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("ul_add_layer", request, adinfo)
                if result is not None:
                    return result
                adata = get_ads().get_adata(adinfo=adinfo)
                layer_name = request.layer_name

                # Check if layer already exists
                if layer_name in adata.layers:
                    raise ValueError(
                        f"Layer '{layer_name}' already exists in adata.layers"
                    )
                # Add the data as a new layer
                adata.layers[layer_name] = adata.X.copy()

                func_kwargs = {"layer_name": layer_name}
                add_op_log(adata, "add_layer", func_kwargs, adinfo)

                return {
                    "status": "success",
                    "message": f"Successfully added layer '{layer_name}' to adata.layers",
                }
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _add_layer, name="add_layer", enabled=True, tags=["preset"]
        )

    def _tool_check_samples(self):
        def _check_samples(
            request: Union[None, str, dict] = None,
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """check the stored samples"""
            # Deserialize parameters - request can be None for this function
            if request is not None:
                request = deserialize_mcp_param(request, type(None))
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                ads = get_ads()
                return {
                    "sampleid": [
                        list(ads.adata_dic[dk].keys()) for dk in ads.adata_dic.keys()
                    ]
                }
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _check_samples, name="check_samples", enabled=True, tags=["preset"]
        )

    def _tool_map_cell_type(self):
        def _map_cell_type(
            request: Union[CelltypeMapCellTypeParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """
            Map cluster IDs to cell type names in the AnnData object.

            This function assigns biological cell type names to clusters identified in adata.obs[cluster_key].
            The mapping can be provided either as a dictionary (mapping cluster IDs to cell type names) or as a list of new category names.
            The result is stored in a new column (added_key) in adata.obs, or the categories of cluster_key are renamed.

            Args:
                request: CelltypeMapCellTypeModel or compatible dict/str, specifying cluster_key, added_key, and mapping or new_names.
                adinfo: AdataInfo or compatible dict/str, identifying the AnnData object.

            Returns:
                list: A list containing a dict with status, message, and updated AnnData object info.
            """
            try:
                # Deserialize parameters
                request = deserialize_mcp_param(request, CelltypeMapCellTypeParam)
                adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())

                result = forward_request("ul_map_cell_type", request, adinfo)
                if result is not None:
                    return result

                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo).copy()
                cluster_key = request.cluster_key
                added_key = request.added_key

                if cluster_key not in adata.obs.columns:
                    raise ValueError(
                        f"cluster key '{cluster_key}' not found in adata.obs"
                    )

                # Perform mapping
                if request.mapping is not None:
                    # Map using provided dictionary
                    adata.obs[added_key] = pd.Categorical(
                        pd.Series(adata.obs[cluster_key]).map(request.mapping)
                    )

                func_kwargs = {
                    "cluster_key": cluster_key,
                    "added_key": added_key,
                    "mapping": request.mapping,
                }
                add_op_log(adata, "map_cell_type", func_kwargs, adinfo)

                # Update the AnnData in the context
                ads.set_adata(adata, adinfo=adinfo)

                return [
                    {
                        "status": "success",
                        "message": f"Successfully mapped values from '{cluster_key}' to '{added_key}'",
                        "adata": adata,
                    }
                ]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _map_cell_type, name="map_cell_type", enabled=True, tags=["preset"]
        )

    def _tool_describe_obs_column(self):
        def _describe_obs_column(
            request: Union[DescribeObsColumnParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """
            Describe a column in adata.obs with appropriate statistics based on data type.

            For numeric (float) columns: returns .obs[name].describe()
            For categorical/string columns with fewer than 500 unique values: returns .obs[name].value_counts()

            Args:
                request: DescribeObsColumnParam or compatible dict/str, specifying the column_name to describe.
                adinfo: AdataInfo or compatible dict/str, identifying the AnnData object.

            Returns:
                dict: Descriptive statistics or value counts for the specified column.
            """
            try:
                # Deserialize parameters
                request = deserialize_mcp_param(request, DescribeObsColumnParam)
                adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())

                result = forward_request("ul_describe_obs_column", request, adinfo)
                if result is not None:
                    return result

                adata = get_ads().get_adata(adinfo=adinfo)
                column_name = request.column_name

                if column_name not in adata.obs.columns:
                    raise ValueError(f"Column '{column_name}' not found in adata.obs")

                column_data = adata.obs[column_name]

                # Check if column is numeric (float)
                if pd.api.types.is_numeric_dtype(column_data):
                    result = column_data.describe().to_dict()
                    result_type = "numeric_describe"
                else:
                    # For categorical or string columns
                    unique_count = column_data.nunique()
                    if unique_count < 500:
                        result = column_data.value_counts().to_dict()
                        result_type = "categorical_value_counts"
                    else:
                        result = {
                            "unique_count": unique_count,
                            "message": f"Column has {unique_count} unique values (>= 500), showing only count",
                        }
                        result_type = "high_cardinality_summary"

                func_kwargs = {"column_name": column_name}
                add_op_log(adata, "describe_obs_column", func_kwargs, adinfo)

                return {
                    "column_name": column_name,
                    "result_type": result_type,
                    "data": result,
                }

            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _describe_obs_column,
            name="describe_obs_column",
            enabled=True,
            tags=["preset"],
        )

    def _tool_run_deseq2_differential_expression(self):
        def _run_deseq2_differential_expression(
            request: Union[DeSeq2DEParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """
            Run pyDeSeq2 differential expression analysis on AnnData object.
            
            This function performs differential gene expression analysis using pyDeSeq2 
            and returns a table of statistically significant differentially expressed genes.
            
            Args:
                request: DeSeq2DEParam or compatible dict/str, specifying the analysis parameters
                adinfo: AdataInfo or compatible dict/str, identifying the AnnData object
                
            Returns:
                pandas.DataFrame: Table of significantly differentially expressed genes with columns:
                    - gene_name: Gene identifier
                    - log2FoldChange: Log2 fold change (group1 vs group2)
                    - pvalue: Raw p-value
                    - padj: Benjamini-Hochberg adjusted p-value
                    - contrast: Contrast description (group1_vs_group2)
                    - condition_column: Name of the condition column used
                    
            Raises:
                ImportError: If pydeseq2 or scanpy is not installed
                ValueError: If input validation fails (missing columns, groups, negative values)
            """
            try:
                # Deserialize parameters
                request = deserialize_mcp_param(request, DeSeq2DEParam)
                adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())

                result = forward_request("ul_run_deseq2_differential_expression", request, adinfo)
                if result is not None:
                    return result

                adata = get_ads().get_adata(adinfo=adinfo)
                
                # Run the analysis
                result_df = self._run_deseq2_analysis(
                    adata=adata,
                    contrast=request.contrast,
                    alpha=request.alpha,
                    lfc_threshold=request.lfc_threshold,
                    n_cpus=request.n_cpus,
                    layer=request.layer
                )
                
                # Log the operation
                func_kwargs = {
                    "contrast": request.contrast,
                    "alpha": request.alpha,
                    "lfc_threshold": request.lfc_threshold,
                    "n_cpus": request.n_cpus,
                    "layer": request.layer
                }
                add_op_log(adata, "run_deseq2_differential_expression", func_kwargs, adinfo)
                
                return result_df
                
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _run_deseq2_differential_expression,
            name="run_deseq2_differential_expression",
            enabled=True,
            tags=["deseq2", "differential_expression", "preset"],
        )

    def _run_deseq2_analysis(
        self,
        adata,
        contrast,
        alpha=0.05,
        lfc_threshold=0.0,
        n_cpus=1,
        layer=None
    ):
        """
        Internal method to run pyDeSeq2 differential expression analysis.
        
        Args:
            adata: AnnData object with count data
            contrast: List of [column_name, group1, group2] for comparison
            alpha: Significance threshold for adjusted p-values (default: 0.05)
            lfc_threshold: Log fold change threshold (default: 0.0)
            n_cpus: Number of CPUs to use (default: 1)
            layer: Layer of AnnData to use for count data. If None, uses adata.X (default: None)
            
        Returns:
            pandas.DataFrame: Results with significantly differentially expressed genes
        """
        import pandas as pd
        import numpy as np
        
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
        except ImportError:
            raise ImportError("pydeseq2 is required for differential expression analysis. Install with: pip install pydeseq2")
        
        try:
            from scanpy.get import _get_obs_rep
        except ImportError:
            raise ImportError("scanpy is required for layer access. Install with: pip install scanpy")
        
        # Validate inputs
        if len(contrast) != 3:
            raise ValueError("contrast must have exactly 3 elements: [column_name, group1, group2]")
        
        condition_col, group1, group2 = contrast
        
        # Check if condition column exists
        if condition_col not in adata.obs.columns:
            raise ValueError(f"Condition column '{condition_col}' not found in adata.obs")
        
        # Check if groups exist in the condition column
        condition_values = adata.obs[condition_col].unique()
        if group1 not in condition_values:
            raise ValueError(f"Group '{group1}' not found in condition column '{condition_col}'. Available groups: {list(condition_values)}")
        if group2 not in condition_values:
            raise ValueError(f"Group '{group2}' not found in condition column '{condition_col}'. Available groups: {list(condition_values)}")
        
        # Get count data from specified layer or X
        count_matrix = _get_obs_rep(adata, layer=layer)
        
        # Check for non-negative count data
        if hasattr(count_matrix, 'min'):
            min_val = count_matrix.min()
        else:
            min_val = np.min(count_matrix)
        
        if min_val < 0:
            layer_desc = f"layer '{layer}'" if layer is not None else "adata.X"
            raise ValueError(f"pyDeSeq2 requires non-negative count data. Found negative values in {layer_desc}")
        
        # Ensure integer counts
        if hasattr(count_matrix, 'astype'):
            count_matrix = count_matrix.astype(int)
        else:
            count_matrix = np.array(count_matrix, dtype=int)
        
        # Create counts DataFrame with genes as columns and samples as rows
        counts_df = pd.DataFrame(
            count_matrix,
            index=adata.obs_names,
            columns=adata.var_names
        )
        
        # Create design matrix
        design_df = adata.obs[[condition_col]].copy()
        
        # Create DeseqDataSet
        try:
            dds = DeseqDataSet(
                counts=counts_df,
                metadata=design_df,
                design=f"~ {condition_col}",
                n_cpus=n_cpus
            )
            
            # Run DESeq2 analysis
            dds.deseq2()
            
            # Get statistics for the contrast
            stat_res = DeseqStats(dds, contrast=[condition_col, group1, group2], n_cpus=n_cpus)
            stat_res.summary()
            
            # Get results for specific contrast
            results_df = stat_res.results_df.copy()
            
            # Filter for significant results
            significant_mask = (
                (results_df['padj'] < alpha) & 
                (results_df['padj'].notna()) &
                (abs(results_df['log2FoldChange']) >= lfc_threshold)
            )
            
            significant_results = results_df[significant_mask].copy()
            
            # Add gene names and contrast information
            significant_results['gene_name'] = significant_results.index
            significant_results['contrast'] = f"{group1}_vs_{group2}"
            significant_results['condition_column'] = condition_col
            
            # Reset index to make gene_name a regular column
            significant_results = significant_results.reset_index(drop=True)
            
            # Sort by adjusted p-value
            significant_results = significant_results.sort_values('padj')
            
            # Select and reorder columns
            output_columns = [
                'gene_name', 'log2FoldChange', 'pvalue', 'padj', 
                'contrast', 'condition_column'
            ]
            
            # Ensure all expected columns exist
            for col in output_columns:
                if col not in significant_results.columns:
                    if col in ['contrast', 'condition_column']:
                        continue  # These are already added above
                    else:
                        significant_results[col] = np.nan
            
            return significant_results[output_columns]
            
        except Exception as e:
            # Handle pyDeSeq2-specific errors gracefully
            if "too small" in str(e).lower() or "insufficient" in str(e).lower():
                # Return empty DataFrame with correct structure for small datasets
                empty_df = pd.DataFrame(columns=[
                    'gene_name', 'log2FoldChange', 'pvalue', 'padj', 
                    'contrast', 'condition_column'
                ])
                return empty_df
            else:
                raise e
