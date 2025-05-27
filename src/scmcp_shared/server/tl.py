from fastmcp import FastMCP, Context
import os
import inspect
from fastmcp import FastMCP
import scanpy as sc
from fastmcp.exceptions import ToolError
from ..schema.tl import *
from ..schema import AdataModel
from scmcp_shared.util import filter_args, add_op_log, forward_request, get_ads, generate_msg


tl_mcp = FastMCP("ScanpyMCP-TL-Server")


@tl_mcp.tool()
async def tsne(
    request: TSNEModel = TSNEModel(),
    adinfo: AdataModel = AdataModel()
):
    """t-distributed stochastic neighborhood embedding (t-SNE) for visualization"""

    try:
        result = await forward_request("tl_tsne", request, adinfo)
        if result is not None:
            return result
        func_kwargs = filter_args(request, sc.tl.tsne)
        ads = get_ads()
        adata = ads.get_adata(adinfo=adinfo)
        sc.tl.tsne(adata, **func_kwargs)
        add_op_log(adata, sc.tl.tsne, func_kwargs, adinfo)
        return generate_msg(adinfo, adata, ads)
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)


@tl_mcp.tool()
async def umap(
    request: UMAPModel = UMAPModel(),
    adinfo: AdataModel = AdataModel()
):
    """Uniform Manifold Approximation and Projection (UMAP) for visualization"""

    try:
        result = await forward_request("tl_umap", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)

@tl_mcp.tool()
async def draw_graph(
    request: DrawGraphModel = DrawGraphModel(),
    adinfo: AdataModel = AdataModel()
):
    """Force-directed graph drawing"""

    try:
        result = await forward_request("tl_draw_graph", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)

@tl_mcp.tool()
async def diffmap(
    request: DiffMapModel = DiffMapModel(),
    adinfo: AdataModel = AdataModel()
):
    """Diffusion Maps for dimensionality reduction"""

    try:
        result = await forward_request("tl_diffmap", request, adinfo)
        if result is not None:
            return result    
        func_kwargs = filter_args(request, sc.tl.diffmap)
        ads = get_ads()
        adata = ads.get_adata(adinfo=adinfo)
        sc.tl.diffmap(adata, **func_kwargs)
        adata.obsm["X_diffmap"] = adata.obsm["X_diffmap"][:,1:]
        add_op_log(adata, sc.tl.diffmap, func_kwargs, adinfo)
        return [generate_msg(adinfo, adata, ads)]
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)

@tl_mcp.tool()
async def embedding_density(
    request: EmbeddingDensityModel = EmbeddingDensityModel(),
    adinfo: AdataModel = AdataModel()
):
    """Calculate the density of cells in an embedding"""

    try:
        result = await forward_request("tl_embedding_density", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)


@tl_mcp.tool()
async def leiden(
    request: LeidenModel = LeidenModel(),
    adinfo: AdataModel = AdataModel()
):
    """Leiden clustering algorithm for community detection"""

    try:
        result = await forward_request("tl_leiden", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)


@tl_mcp.tool()
async def louvain(
    request: LouvainModel = LouvainModel(),
    adinfo: AdataModel = AdataModel()
):
    """Louvain clustering algorithm for community detection"""

    try:
        result = await forward_request("tl_louvain", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)

@tl_mcp.tool()
async def dendrogram(
    request: DendrogramModel,
    adinfo: AdataModel = AdataModel()
):
    """Hierarchical clustering dendrogram"""

    try:
        result = await forward_request("tl_dendrogram", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)

@tl_mcp.tool()
async def dpt(
    request: DPTModel = DPTModel(),
    adinfo: AdataModel = AdataModel()
):
    """Diffusion Pseudotime (DPT) analysis"""

    try:
        result = await forward_request("tl_dpt", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)


@tl_mcp.tool()
async def paga(
    request: PAGAModel = PAGAModel(),
    adinfo: AdataModel = AdataModel()
):
    """Partition-based graph abstraction"""

    try:
        result = await forward_request("tl_paga", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)


@tl_mcp.tool()
async def ingest(
    request: IngestModel = IngestModel(),
    adinfo: AdataModel = AdataModel()
):
    """Map labels and embeddings from reference data to new data"""

    try:
        result = await forward_request("tl_ingest", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)

@tl_mcp.tool()
async def rank_genes_groups(
    request: RankGenesGroupsModel,
    adinfo: AdataModel = AdataModel()
):
    """Rank genes for characterizing groups, for differentially expressison analysis"""

    try:
        result = await forward_request("tl_rank_genes_groups", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)


@tl_mcp.tool()
async def filter_rank_genes_groups(
    request: FilterRankGenesGroupsModel = FilterRankGenesGroupsModel(),
    adinfo: AdataModel = AdataModel()
):
    """Filter out genes based on fold change and fraction of genes"""

    try:
        result = await forward_request("tl_filter_rank_genes_groups", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)


@tl_mcp.tool()
async def marker_gene_overlap(
    request: MarkerGeneOverlapModel = MarkerGeneOverlapModel(),
    adinfo: AdataModel = AdataModel()
):
    """Calculate overlap between data-derived marker genes and reference markers"""

    try:
        result = await forward_request("tl_marker_gene_overlap", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)

@tl_mcp.tool()
async def score_genes(
    request: ScoreGenesModel,
    adinfo: AdataModel = AdataModel()
):
    """Score a set of genes based on their average expression"""
    try:
        result = await forward_request("tl_score_genes", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)

@tl_mcp.tool()
async def score_genes_cell_cycle(
    request: ScoreGenesCellCycleModel,
    adinfo: AdataModel = AdataModel()
):
    """Score cell cycle genes and assign cell cycle phases"""

    try:
        result = await forward_request("tl_score_genes_cell_cycle", request, adinfo)
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
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)


@tl_mcp.tool()
async def pca(
    request: PCAModel = PCAModel(),
    adinfo: AdataModel = AdataModel()
):
    """Principal component analysis"""

    try:
        result = await forward_request("tl_pca", request, adinfo)
        if result is not None:
            return result
        func_kwargs = filter_args(request, sc.pp.pca)
        ads = get_ads()
        adata = ads.get_adata(adinfo=adinfo)
        sc.pp.pca(adata, **func_kwargs)
        add_op_log(adata, sc.pp.pca, func_kwargs, adinfo)
        return [
            generate_msg(adinfo, adata, ads)
        ]
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)



mcp1 = FastMCP("ScanpyMCP-TL-Server")

class ScanpyToolsMCP:
    def __init__(self, include_tools: list = None, exclude_tools: list = None):
        """
        Initialize ScanpyMCP with optional tool filtering.
        
        Args:
            include_tools (list, optional): List of tool names to include. If None, all tools are included.
            exclude_tools (list, optional): List of tool names to exclude. If None, no tools are excluded.
        """
        self.mcp = mcp1
        self.include_tools = include_tools
        self.exclude_tools = exclude_tools
        self._register_tools()

    def _register_tools(self):
        """Register all tool methods with the FastMCP instance based on include/exclude filters"""
        # Get all methods of the class
        methods = inspect.getmembers(self, predicate=inspect.ismethod)
        
        # Filter methods that start with _tool_
        tool_methods = {
            name[6:]: method  # Remove '_tool_' prefix
            for name, method in methods
            if name.startswith('_tool_')
        }
        
        # Filter tools based on include/exclude lists
        if self.include_tools is not None:
            tool_methods = {k: v for k, v in tool_methods.items() if k in self.include_tools}
        
        if self.exclude_tools is not None:
            tool_methods = {k: v for k, v in tool_methods.items() if k not in self.exclude_tools}

        # Register filtered tools
        for tool_name, tool_func in tool_methods.items():
            self.mcp.add_tool(tool_func,name=tool_name)

    async def _tool_tsne(
        self,
        request: TSNEModel = TSNEModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """t-distributed stochastic neighborhood embedding (t-SNE) for visualization"""
        try:
            result = await forward_request("tl_tsne", request, adinfo)
            if result is not None:
                return result
            func_kwargs = filter_args(request, sc.tl.tsne)
            ads = get_ads()
            adata = ads.get_adata(adinfo=adinfo)
            sc.tl.tsne(adata, **func_kwargs)
            add_op_log(adata, sc.tl.tsne, func_kwargs, adinfo)
            return generate_msg(adinfo, adata, ads)
        except ToolError as e:
            raise ToolError(e)
        except Exception as e:
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_umap(
        self,
        request: UMAPModel = UMAPModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Uniform Manifold Approximation and Projection (UMAP) for visualization"""
        try:
            result = await forward_request("tl_umap", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_draw_graph(
        self,
        request: DrawGraphModel = DrawGraphModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Force-directed graph drawing"""
        try:
            result = await forward_request("tl_draw_graph", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_diffmap(
        self,
        request: DiffMapModel = DiffMapModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Diffusion Maps for dimensionality reduction"""
        try:
            result = await forward_request("tl_diffmap", request, adinfo)
            if result is not None:
                return result    
            func_kwargs = filter_args(request, sc.tl.diffmap)
            ads = get_ads()
            adata = ads.get_adata(adinfo=adinfo)
            sc.tl.diffmap(adata, **func_kwargs)
            adata.obsm["X_diffmap"] = adata.obsm["X_diffmap"][:,1:]
            add_op_log(adata, sc.tl.diffmap, func_kwargs, adinfo)
            return [generate_msg(adinfo, adata, ads)]
        except ToolError as e:
            raise ToolError(e)
        except Exception as e:
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_embedding_density(
        self,
        request: EmbeddingDensityModel = EmbeddingDensityModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Calculate the density of cells in an embedding"""
        try:
            result = await forward_request("tl_embedding_density", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_leiden(
        self,
        request: LeidenModel = LeidenModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Leiden clustering algorithm for community detection"""
        try:
            result = await forward_request("tl_leiden", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_louvain(
        self,
        request: LouvainModel = LouvainModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Louvain clustering algorithm for community detection"""
        try:
            result = await forward_request("tl_louvain", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_dendrogram(
        self,
        request: DendrogramModel,
        adinfo: AdataModel = AdataModel()
    ):
        """Hierarchical clustering dendrogram"""
        try:
            result = await forward_request("tl_dendrogram", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_dpt(
        self,
        request: DPTModel = DPTModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Diffusion Pseudotime (DPT) analysis"""
        try:
            result = await forward_request("tl_dpt", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_paga(
        self,
        request: PAGAModel = PAGAModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Partition-based graph abstraction"""
        try:
            result = await forward_request("tl_paga", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_ingest(
        self,
        request: IngestModel = IngestModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Map labels and embeddings from reference data to new data"""
        try:
            result = await forward_request("tl_ingest", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_rank_genes_groups(
        self,
        request: RankGenesGroupsModel,
        adinfo: AdataModel = AdataModel()
    ):
        """Rank genes for characterizing groups, for differentially expressison analysis"""
        try:
            result = await forward_request("tl_rank_genes_groups", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_filter_rank_genes_groups(
        self,
        request: FilterRankGenesGroupsModel = FilterRankGenesGroupsModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Filter out genes based on fold change and fraction of genes"""
        try:
            result = await forward_request("tl_filter_rank_genes_groups", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_marker_gene_overlap(
        self,
        request: MarkerGeneOverlapModel = MarkerGeneOverlapModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Calculate overlap between data-derived marker genes and reference markers"""
        try:
            result = await forward_request("tl_marker_gene_overlap", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_score_genes(
        self,
        request: ScoreGenesModel,
        adinfo: AdataModel = AdataModel()
    ):
        """Score a set of genes based on their average expression"""
        try:
            result = await forward_request("tl_score_genes", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_score_genes_cell_cycle(
        self,
        request: ScoreGenesCellCycleModel,
        adinfo: AdataModel = AdataModel()
    ):
        """Score cell cycle genes and assign cell cycle phases"""
        try:
            result = await forward_request("tl_score_genes_cell_cycle", request, adinfo)
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
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)

    async def _tool_pca(
        self,
        request: PCAModel = PCAModel(),
        adinfo: AdataModel = AdataModel()
    ):
        """Principal component analysis"""
        try:
            result = await forward_request("tl_pca", request, adinfo)
            if result is not None:
                return result
            func_kwargs = filter_args(request, sc.pp.pca)
            ads = get_ads()
            adata = ads.get_adata(adinfo=adinfo)
            sc.pp.pca(adata, **func_kwargs)
            add_op_log(adata, sc.pp.pca, func_kwargs, adinfo)
            return [
                generate_msg(adinfo, adata, ads)
            ]
        except ToolError as e:
            raise ToolError(e)
        except Exception as e:
            if hasattr(e, '__context__') and e.__context__:
                raise ToolError(e.__context__)
            else:
                raise ToolError(e)


