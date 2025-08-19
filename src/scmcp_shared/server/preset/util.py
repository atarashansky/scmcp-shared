from typing import Union
from fastmcp.exceptions import ToolError
from fastmcp.tools.tool import Tool
from scmcp_shared.schema.preset.util import *
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
            print(request)
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
