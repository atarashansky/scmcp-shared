import scanpy as sc
import decoupler as dc
from typing import Union
from fastmcp.exceptions import ToolError
from fastmcp.tools.tool import Tool
from scmcp_shared.schema.preset.tl import *
from scmcp_shared.schema.preset import AdataInfo
from scmcp_shared.util import (
    filter_args,
    add_op_log,
    forward_request,
    get_ads,
    generate_msg,
    deserialize_mcp_param,
)
from scmcp_shared.mcp_base import BaseMCP


class ScanpyToolsMCP(BaseMCP):
    def __init__(
        self,
        include_tools: list = None,
        exclude_tools: list = None,
        AdataInfo: AdataInfo = AdataInfo,
    ):
        """
        Initialize ScanpyMCP with optional tool filtering.

        Args:
            include_tools (list, optional): List of tool names to include. If None, all tools are included.
            exclude_tools (list, optional): List of tool names to exclude. If None, no tools are excluded.
            AdataInfo: The AdataInfo class to use for type annotations.
        """
        super().__init__("ScanpyMCP-TL-Server", include_tools, exclude_tools, AdataInfo)

    def _tool_pathway_activity(self):
        def _pathway_activity(
            request: Union[PathwayActivityParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """
            Perform pathway activity inference on single-cell data using decoupler's Multivariate Linear Model (MLM) method.

            This tool leverages the PROGENy gene sets to estimate pathway activities for each sample in the provided AnnData object.
            The MLM method fits a multivariate linear model where observed molecular readouts are the response variable and pathway
            weights are the covariates. The resulting t-values represent the activity of each pathway.

            Parameters:
                request (PathwayActivityParam | str | dict): Parameters for pathway activity inference, including organism, top genes, and MLM options.
                adinfo (AdataInfo | str | dict, optional): Information to locate and identify the AnnData object.

            Returns:
                list[dict]: A list containing dictionaries for the original and pathway activity AnnData objects, with associated sample IDs and types.

            Raises:
                ToolError: If an error occurs during pathway activity inference.
            """
            # Deserialize parameters
            request = deserialize_mcp_param(request, PathwayActivityParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("if_pathway_activity", request, adinfo)
                if result is not None:
                    return result
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                kwargs = (
                    request.model_dump()
                    if hasattr(request, "model_dump")
                    else dict(request)
                )
                progeny = dc.op.progeny(
                    organism=kwargs["organism"], top=kwargs.get("top", None)
                )
                func_kwargs = filter_args(request, dc.mt.mlm)
                dc.mt.mlm(data=adata, net=progeny, **func_kwargs)
                score = dc.pp.get_obsm(adata=adata, key="score_mlm")
                add_op_log(adata, dc.mt.mlm, func_kwargs, adinfo)
                sdtype = "activity"
                ads.set_adata(score, sampleid="score_mlm", sdtype=sdtype)
                return [
                    {
                        "sampleid": getattr(adinfo, "sampleid", None) or ads.active_id,
                        "adtype": getattr(adinfo, "adtype", None),
                        "adata": adata,
                    },
                    {"sampleid": "score_mlm", "adtype": sdtype, "adata": score},
                ]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _pathway_activity,
            name="pathway_activity",
            enabled=True,
            tags=["preset"],
        )

    def _tool_tf_activity(self):
        def _tf_activity(
            request: Union[TFActivityParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """
            Infers transcription factor (TF) activity using decoupler's Univariate Linear Model (ULM) method.

            This tool estimates the activity of transcription factors for each sample by fitting a linear model
            where observed molecular readouts are the response variable and regulator weights are the explanatory variable.
            Target features with no associated weight are set to zero. The resulting t-value from the fitted model
            represents the activity of a given regulator (TF).

            Parameters
            ----------
            request : TFActivityParam, str, or dict
                Parameters for the ULM method, including organism, whether to remove complexes, and license type.
            adinfo : AdataInfo, str, or dict, optional
                Information about the AnnData object to use.

            Returns
            -------
            list of dict
                Contains the original and the TF activity-inferred AnnData objects, with associated sample and adata types.
            """
            # Deserialize parameters
            request = deserialize_mcp_param(request, TFActivityParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("if_tf_activity", request, adinfo)
                if result is not None:
                    return result
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                kwargs = (
                    request.model_dump()
                    if hasattr(request, "model_dump")
                    else dict(request)
                )
                net = dc.op.collectri(organism=kwargs["organism"])
                func_kwargs = filter_args(request, dc.mt.ulm)
                dc.mt.ulm(data=adata, net=net, **func_kwargs)
                score = dc.pp.get_obsm(adata=adata, key="score_ulm")
                add_op_log(adata, dc.mt.ulm, func_kwargs, adinfo)
                sdtype = "activity"
                ads.set_adata(score, sampleid="score_ulm", sdtype=sdtype)
                return [
                    {
                        "sampleid": getattr(adinfo, "sampleid", None) or ads.active_id,
                        "adtype": getattr(adinfo, "adtype", None),
                        "adata": adata,
                    },
                    {"sampleid": "score_ulm", "adtype": sdtype, "adata": score},
                ]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _tf_activity,
            name="tf_activity",
            enabled=True,
            tags=["preset"],
        )

    def _tool_tsne(self):
        def _tsne(
            request: Union[TSNEParam, str, dict] = None,
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """t-distributed stochastic neighborhood embedding (t-SNE) for visualization"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, TSNEParam, TSNEParam())
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_tsne", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.tsne)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.tsne(adata, **func_kwargs)
                add_op_log(adata, sc.tl.tsne, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(_tsne, name="tsne", enabled=True, tags=["preset"])

    def _tool_umap(self):
        def _umap(
            request: Union[UMAPParam, str, dict] = None,
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Uniform Manifold Approximation and Projection (UMAP) for visualization"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, UMAPParam, UMAPParam())
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_umap", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.umap)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.umap(adata, **func_kwargs)
                add_op_log(adata, sc.tl.umap, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(_umap, name="umap", enabled=True, tags=["preset"])

    def _tool_draw_graph(self):
        def _draw_graph(
            request: Union[DrawGraphParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Force-directed graph drawing"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, DrawGraphParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_draw_graph", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.draw_graph)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.draw_graph(adata, **func_kwargs)
                add_op_log(adata, sc.tl.draw_graph, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _draw_graph, name="draw_graph", enabled=True, tags=["preset"]
        )

    def _tool_diffmap(self):
        def _diffmap(
            request: Union[DiffMapParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Diffusion Maps for dimensionality reduction"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, DiffMapParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_diffmap", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.diffmap)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.diffmap(adata, **func_kwargs)
                adata.obsm["X_diffmap"] = adata.obsm["X_diffmap"][:, 1:]
                add_op_log(adata, sc.tl.diffmap, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _diffmap, name="diffmap", enabled=True, tags=["preset"]
        )

    def _tool_embedding_density(self):
        def _embedding_density(
            request: Union[EmbeddingDensityParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Calculate the density of cells in an embedding"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, EmbeddingDensityParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_embedding_density", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.embedding_density)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.embedding_density(adata, **func_kwargs)
                add_op_log(adata, sc.tl.embedding_density, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _embedding_density, name="embedding_density", enabled=True, tags=["preset"]
        )

    def _tool_leiden(self):
        def _leiden(
            request: Union[LeidenParam, str, dict] = None,
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Leiden clustering algorithm for community detection"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, LeidenParam, LeidenParam())
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_leiden", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.leiden)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.leiden(adata, **func_kwargs)
                add_op_log(adata, sc.tl.leiden, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(_leiden, name="leiden", enabled=True, tags=["preset"])

    def _tool_louvain(self):
        def _louvain(
            request: Union[LouvainParam, str, dict] = None,
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Louvain clustering algorithm for community detection"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, LouvainParam, LouvainParam())
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_louvain", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.louvain)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.louvain(adata, **func_kwargs)
                add_op_log(adata, sc.tl.louvain, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _louvain, name="louvain", enabled=True, tags=["preset"]
        )

    def _tool_dendrogram(self):
        def _dendrogram(
            request: Union[DendrogramParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Hierarchical clustering dendrogram"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, DendrogramParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_dendrogram", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.dendrogram)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.dendrogram(adata, **func_kwargs)
                add_op_log(adata, sc.tl.dendrogram, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _dendrogram, name="dendrogram", enabled=True, tags=["preset"]
        )

    def _tool_dpt(self):
        def _dpt(
            request: Union[DPTParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Diffusion Pseudotime (DPT) analysis"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, DPTParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_dpt", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.dpt)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.dpt(adata, **func_kwargs)
                add_op_log(adata, sc.tl.dpt, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(_dpt, name="dpt", enabled=True, tags=["preset"])

    def _tool_paga(self):
        def _paga(
            request: Union[PAGAParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Partition-based graph abstraction"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, PAGAParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_paga", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.paga)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.paga(adata, **func_kwargs)
                add_op_log(adata, sc.tl.paga, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(_paga, name="paga", enabled=True, tags=["preset"])

    def _tool_ingest(self):
        def _ingest(
            request: Union[IngestParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Map labels and embeddings from reference data to new data"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, IngestParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_ingest", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.ingest)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.ingest(adata, **func_kwargs)
                add_op_log(adata, sc.tl.ingest, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(_ingest, name="ingest", enabled=True, tags=["preset"])

    def _tool_rank_genes_groups(self):
        def _rank_genes_groups(
            request: Union[RankGenesGroupsParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Rank genes for characterizing groups, for differentially expressison analysis"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, RankGenesGroupsParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_rank_genes_groups", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.rank_genes_groups)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.rank_genes_groups(adata, **func_kwargs)
                add_op_log(adata, sc.tl.rank_genes_groups, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _rank_genes_groups, name="rank_genes_groups", enabled=True, tags=["preset"]
        )

    def _tool_filter_rank_genes_groups(self):
        def _filter_rank_genes_groups(
            request: Union[FilterRankGenesGroupsParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Filter out genes based on fold change and fraction of genes"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, FilterRankGenesGroupsParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_filter_rank_genes_groups", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.filter_rank_genes_groups)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.filter_rank_genes_groups(adata, **func_kwargs)
                add_op_log(adata, sc.tl.filter_rank_genes_groups, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _filter_rank_genes_groups,
            name="filter_rank_genes_groups",
            enabled=True,
            tags=["preset"],
        )

    def _tool_marker_gene_overlap(self):
        def _marker_gene_overlap(
            request: Union[MarkerGeneOverlapParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Calculate overlap between data-derived marker genes and reference markers"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, MarkerGeneOverlapParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_marker_gene_overlap", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.marker_gene_overlap)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.marker_gene_overlap(adata, **func_kwargs)
                add_op_log(adata, sc.tl.marker_gene_overlap, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _marker_gene_overlap,
            name="marker_gene_overlap",
            enabled=True,
            tags=["preset"],
        )

    def _tool_score_genes(self):
        def _score_genes(
            request: Union[ScoreGenesParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Score a set of genes based on their average expression"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, ScoreGenesParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_score_genes", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.score_genes)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.score_genes(adata, **func_kwargs)
                add_op_log(adata, sc.tl.score_genes, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _score_genes, name="score_genes", enabled=True, tags=["preset"]
        )

    def _tool_score_genes_cell_cycle(self):
        def _score_genes_cell_cycle(
            request: Union[ScoreGenesCellCycleParam, str, dict],
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Score cell cycle genes and assign cell cycle phases"""
            # Deserialize parameters
            request = deserialize_mcp_param(request, ScoreGenesCellCycleParam)
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_score_genes_cell_cycle", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.score_genes_cell_cycle)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.score_genes_cell_cycle(adata, **func_kwargs)
                add_op_log(adata, sc.tl.score_genes_cell_cycle, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(
            _score_genes_cell_cycle,
            name="score_genes_cell_cycle",
            enabled=True,
            tags=["preset"],
        )

    def _tool_pca(self):
        def _pca(
            request: Union[PCAParam, str, dict] = None,
            adinfo: Union[AdataInfo, str, dict] = None,
        ):
            """Compute PCA (Principal Component Analysis)."""
            # Deserialize parameters
            request = deserialize_mcp_param(request, PCAParam, PCAParam())
            adinfo = deserialize_mcp_param(adinfo, self.AdataInfo, self.AdataInfo())
            try:
                result = forward_request("tl_pca", request, adinfo)
                if result is not None:
                    return result
                func_kwargs = filter_args(request, sc.tl.pca)
                ads = get_ads()
                adata = ads.get_adata(adinfo=adinfo)
                sc.tl.pca(adata, **func_kwargs)
                add_op_log(adata, sc.tl.pca, func_kwargs, adinfo)
                return [generate_msg(adinfo, adata, ads)]
            except ToolError as e:
                raise ToolError(e)
            except Exception as e:
                if hasattr(e, "__context__") and e.__context__:
                    raise ToolError(e.__context__)
                else:
                    raise ToolError(e)

        return Tool.from_function(_pca, name="pca", enabled=True, tags=["preset"])
