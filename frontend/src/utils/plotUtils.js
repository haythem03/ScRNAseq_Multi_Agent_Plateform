export function collectBase64Plots(payload = {}, stepId = 'qc') {
    const candidates = [
        payload?.plots?.[stepId],
        payload?.steps?.[stepId]?.plots,
        payload?.steps?.[stepId]?.plot,
        payload?.steps?.[stepId]?.result?.plots,
        payload?.steps?.[stepId]?.result?.plot,
        payload?.[stepId]?.plots,
        payload?.[stepId]?.plot,
        payload?.plot,
        payload?.plots,
        payload?.qc_results?.plots,
        payload?.qc_results?.plot,
    ];

    const aggregated = {};

    candidates.forEach((candidate) => {
        if (!candidate) return;
        if (typeof candidate === 'string') {
            const key = Object.keys(aggregated).length === 0 ? 'plot' : `plot_${Object.keys(aggregated).length + 1}`;
            aggregated[key] = candidate;
            return;
        }
        if (typeof candidate === 'object') {
            Object.entries(candidate).forEach(([key, value]) => {
                if (typeof value === 'string' && value.length > 20) {
                    aggregated[key] = aggregated[key] || value;
                }
            });
        }
    });

    return aggregated;
}

export function prettifyPlotLabel(key = '') {
    return key
        .replace(/[_-]+/g, ' ')
        .replace(/\b\w/g, (chr) => chr.toUpperCase())
        .trim() || 'Plot';
}
