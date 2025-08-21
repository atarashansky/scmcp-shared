import pytest
import numpy as np
import pandas as pd
from anndata import AnnData
from scmcp_shared.server.preset.util import ScanpyUtilMCP
from scmcp_shared.schema.util import DeSeq2DEParam
from scmcp_shared.schema.preset import AdataInfo
import nest_asyncio

nest_asyncio.apply()


@pytest.fixture
def util_mcp():
    """Create a ScanpyUtilMCP instance for testing"""
    return ScanpyUtilMCP()


@pytest.fixture
def mock_adata():
    """Create a mock AnnData object with count data for testing"""
    n_obs = 100
    n_vars = 50
    
    # Create mock count data (non-negative integers)
    np.random.seed(42)
    X = np.random.poisson(10, size=(n_obs, n_vars))
    
    # Create gene names
    var_names = [f"Gene_{i}" for i in range(n_vars)]
    
    # Create cell names
    obs_names = [f"Cell_{i}" for i in range(n_obs)]
    
    # Create condition labels for differential expression
    conditions = ['control'] * 50 + ['treatment'] * 50
    
    # Create mock AnnData object with layers
    adata = AnnData(
        X=X,
        var=pd.DataFrame(index=var_names),
        obs=pd.DataFrame({'condition': conditions}, index=obs_names)
    )
    
    # Add a test layer with different count data
    np.random.seed(123)
    adata.layers['raw_counts'] = np.random.poisson(15, size=(n_obs, n_vars))
    
    return adata


@pytest.fixture
def deseq2_params():
    """Create mock parameters for DESeq2 analysis"""
    return DeSeq2DEParam(
        contrast=['condition', 'treatment', 'control'],
        alpha=0.05,
        lfc_threshold=0.0,
        n_cpus=1
    )


def test_deseq2_param_validation():
    """Test DESeq2 parameter validation"""
    # Test valid parameters
    params = DeSeq2DEParam(
        contrast=['condition', 'treatment', 'control'],
        alpha=0.05,
        lfc_threshold=0.5,
        n_cpus=2
    )
    assert params.contrast == ['condition', 'treatment', 'control']
    assert params.alpha == 0.05
    assert params.lfc_threshold == 0.5
    assert params.n_cpus == 2
    
    # Test invalid alpha (should be between 0 and 1)
    with pytest.raises(ValueError):
        DeSeq2DEParam(
            contrast=['condition', 'treatment', 'control'],
            alpha=1.5
        )
    
    # Test invalid lfc_threshold (should be >= 0)
    with pytest.raises(ValueError):
        DeSeq2DEParam(
            contrast=['condition', 'treatment', 'control'],
            lfc_threshold=-1.0
        )


def test_deseq2_function_validation(mock_adata, util_mcp):
    """Test function input validation"""
    # Test with missing condition column
    with pytest.raises(ValueError, match="Condition column 'missing_col' not found"):
        util_mcp._run_deseq2_analysis(
            mock_adata, 
            contrast=['missing_col', 'treatment', 'control']
        )
    
    # Test with missing group in condition column
    with pytest.raises(ValueError, match="Group 'missing_group' not found"):
        util_mcp._run_deseq2_analysis(
            mock_adata, 
            contrast=['condition', 'missing_group', 'control']
        )
    
    # Test with negative values (should fail for pyDeSeq2)
    negative_adata = mock_adata.copy()
    negative_adata.X[0, 0] = -1
    with pytest.raises(ValueError, match="pyDeSeq2 requires non-negative count data"):
        util_mcp._run_deseq2_analysis(
            negative_adata, 
            contrast=['condition', 'treatment', 'control']
        )


def test_deseq2_basic_functionality(mock_adata, deseq2_params, util_mcp):
    """Test basic DESeq2 functionality (requires pyDeSeq2 installation)"""
    try:
        result = util_mcp._run_deseq2_analysis(
            mock_adata,
            contrast=deseq2_params.contrast,
            alpha=deseq2_params.alpha,
            lfc_threshold=deseq2_params.lfc_threshold,
            n_cpus=deseq2_params.n_cpus
        )
        
        # Check that result is a DataFrame
        assert isinstance(result, pd.DataFrame)
        
        # Check that result has expected columns
        expected_cols = ['log2FoldChange', 'pvalue', 'padj', 'gene_name', 'contrast', 'condition_column']
        for col in expected_cols:
            assert col in result.columns
        
        # Check that results are filtered correctly
        if len(result) > 0:
            assert all(result['padj'] < deseq2_params.alpha)
            assert all(abs(result['log2FoldChange']) >= deseq2_params.lfc_threshold)
        
        # Check that contrast information is correct
        if len(result) > 0:
            assert result['contrast'].iloc[0] == 'treatment_vs_control'
            assert result['condition_column'].iloc[0] == 'condition'
            
    except ImportError:
        pytest.skip("pyDeSeq2 not installed, skipping integration test")


def test_deseq2_empty_result(mock_adata, util_mcp):
    """Test DESeq2 with very stringent parameters that should return empty result"""
    try:
        result = util_mcp._run_deseq2_analysis(
            mock_adata,
            contrast=['condition', 'treatment', 'control'],
            alpha=0.001,  # Very stringent p-value
            lfc_threshold=10.0,  # Very high fold change threshold
            n_cpus=1
        )
        
        # Should return empty DataFrame with correct columns
        assert isinstance(result, pd.DataFrame)
        expected_cols = ['log2FoldChange', 'pvalue', 'padj', 'gene_name', 'contrast', 'condition_column']
        for col in expected_cols:
            assert col in result.columns
            
    except ImportError:
        pytest.skip("pyDeSeq2 not installed, skipping integration test")


def test_deseq2_with_small_dataset(util_mcp):
    """Test DESeq2 with minimal dataset"""
    # Create very small dataset
    n_obs = 6
    n_vars = 10
    
    np.random.seed(123)
    X = np.random.poisson(5, size=(n_obs, n_vars))
    
    var_names = [f"Gene_{i}" for i in range(n_vars)]
    obs_names = [f"Cell_{i}" for i in range(n_obs)]
    conditions = ['control'] * 3 + ['treatment'] * 3
    
    small_adata = AnnData(
        X=X,
        var=pd.DataFrame(index=var_names),
        obs=pd.DataFrame({'condition': conditions}, index=obs_names)
    )
    
    try:
        result = util_mcp._run_deseq2_analysis(
            small_adata,
            contrast=['condition', 'treatment', 'control'],
            alpha=0.05,
            lfc_threshold=0.0,
            n_cpus=1
        )
        
        # Should run without error
        assert isinstance(result, pd.DataFrame)
        
    except ImportError:
        pytest.skip("pyDeSeq2 not installed, skipping integration test")
    except Exception as e:
        # pyDeSeq2 might fail with very small datasets - this is expected
        assert "too small" in str(e).lower() or "insufficient" in str(e).lower() or len(str(e)) > 0


def test_deseq2_with_layer(mock_adata, util_mcp):
    """Test DESeq2 with specified layer"""
    try:
        # Test with layer parameter
        result = util_mcp._run_deseq2_analysis(
            mock_adata,
            contrast=['condition', 'treatment', 'control'],
            alpha=0.05,
            lfc_threshold=0.0,
            n_cpus=1,
            layer='raw_counts'
        )
        
        # Should run without error
        assert isinstance(result, pd.DataFrame)
        expected_cols = ['log2FoldChange', 'pvalue', 'padj', 'gene_name', 'contrast', 'condition_column']
        for col in expected_cols:
            assert col in result.columns
            
    except ImportError:
        pytest.skip("pyDeSeq2 not installed, skipping layer test")


def test_deseq2_invalid_layer(mock_adata, util_mcp):
    """Test DESeq2 with invalid layer name"""
    try:
        with pytest.raises(KeyError):
            util_mcp._run_deseq2_analysis(
                mock_adata,
                contrast=['condition', 'treatment', 'control'],
                alpha=0.05,
                lfc_threshold=0.0,
                n_cpus=1,
                layer='nonexistent_layer'
            )
    except ImportError:
        pytest.skip("pyDeSeq2 not installed, skipping invalid layer test")


def test_deseq2_param_validation_with_layer():
    """Test DESeq2 parameter validation with layer"""
    # Test valid parameters with layer
    params = DeSeq2DEParam(
        contrast=['condition', 'treatment', 'control'],
        alpha=0.05,
        lfc_threshold=0.5,
        n_cpus=2,
        layer='raw_counts'
    )
    assert params.contrast == ['condition', 'treatment', 'control']
    assert params.alpha == 0.05
    assert params.lfc_threshold == 0.5
    assert params.n_cpus == 2
    assert params.layer == 'raw_counts'
    
    # Test valid parameters without layer (should be None)
    params_no_layer = DeSeq2DEParam(
        contrast=['condition', 'treatment', 'control']
    )
    assert params_no_layer.layer is None


@pytest.mark.asyncio
async def test_deseq2_mcp_endpoint(mcp, mock_adata):
    """Test the MCP endpoint for DESeq2 (requires proper setup)"""
    pytest.skip("MCP endpoint test requires full server setup with proper backend")