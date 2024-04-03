import pandas as pd
import altair as alt
import datapane as dp

class Report:
    def __init__(self, data, path="report.html", open=False):
        """
        Initialize the Report class with data, output path, and open flag.

        Parameters:
        - data: The data to be loaded into the report.
        - path: The output path for the report file. Default is "report.html".
        - open: Flag indicating whether to open the report file after creation. Default is False.
        """
        self._load_data(data)
        self.path = path
        self.open = open
        self.d = 'd_call' in self.data.columns
        self.non_productive = any(self.data['productive'] == False)
        self.corrupted = any(self.data['v_sequence_start'] > 0)

    def _load_data(self, data):
        if isinstance(data, list):
            self.data = pd.DataFrame(data)
        elif isinstance(data, pd.DataFrame):
            self.data = data
        else:
            raise ValueError("Invalid data type. Expected a list or a DataFrame.")
        
    def create_report(self):
        """
        Create the report layout and save the report.
        """
        self._report_layout()
        self._save_report()
    
    def _report_layout(self):
        """
        Define the layout of the report.
        """
        banner_html = """<div style="padding: 10px;display: flex;align-items: center;font-size: 40px;color: #312E81;background: #EEF2FF;">
        <h1>OverView of simulated data using GenAIRR</h1>
        </div>
        """
        header = dp.Group(dp.HTML(banner_html), columns=1)

        self.layout = dp.View(
            header,
            dp.Group(
                self._nreads(),
                self._percentage_productive(),
                columns=2
            ),
            dp.Text("Distribution of Allele Usage"),
            self._allele_distribution_plot_group(),
            dp.Text("Reads productivity"),
            self._productive_distribution_plot_group(),
            dp.Text("General statistics"),
            self._general_statistics_group()
        )

    def _nreads(self):
        """
        Generate number of reads Big Number component.
        """
        nreads = dp.BigNumber(
            heading="Number of reads",
            value=f"{self.data.shape[0]:,}"
        )
        return nreads

    def _percentage_productive(self):
        """
        Generate percentage of productive reads Big Number component.
        """
        counts = sum(self.data['productive'] == True) / self.data.shape[0] * 100
        counts = str(round(counts, 2)) + '%'
        return dp.BigNumber(
            heading="Percentage of productive reads",
            value=counts
        )
    
    def _allele_distribution_plot(self, call):
        """
        Generate an allele distribution plot for a given call.
        """
        data = self._allele_distribution(call)
        chart = alt.Chart(data).mark_bar().encode(
            x=call,
            y=alt.Y('count', axis=alt.Axis(format='%'), title='Percentage of reads'),
            tooltip=[alt.Tooltip(call, title='Allele Call')]
        )
        return chart

    def _allele_distribution(self, call):
        """
        Compute allele distribution data.
        """
        tmp = self.data[call].apply(lambda x: x.split(',')[0])
        allele_counts = tmp.value_counts()
        total = allele_counts.sum()
        allele_usage = allele_counts / total
        allele_usage = allele_usage.to_frame().reset_index().rename(columns={'index': call, call: 'count'})
        return allele_usage
    
    def _allele_distribution_plot_group(self):
        """
        Generate allele distribution plot group.
        """
        if self.d:
            gr_plot = dp.View(
                dp.Plot(self._allele_distribution_plot('v_call').properties(width=1800)),
                dp.Group(
                    dp.Plot(self._allele_distribution_plot('d_call')),
                    dp.Plot(self._allele_distribution_plot('j_call')),
                    columns=2
                )
            )
        else:
            gr_plot = dp.View(
                dp.Plot(self._allele_distribution_plot('v_call').properties(width=1800)),
                dp.Plot(self._allele_distribution_plot('j_call'))
            )
        return gr_plot

    def _productive_distribution_plot(self):
        """
        Generate a productive distribution plot.
        """
        data = self._productive_distribution()
        chart = alt.Chart(data).mark_bar().encode(
            x=alt.X('productive', axis=alt.Axis(labelAngle=0)),
            y=alt.Y('count', axis=alt.Axis(format='%'), title='Percentage of reads'),
            color=alt.Color('color').scale(None)
        )
        return chart
    
    def _productive_distribution(self):
        """
        Compute productive distribution data.
        """
        productive_counts = self.data['productive'].value_counts()
        total = productive_counts.sum()
        productive_usage = productive_counts / total
        productive_usage = productive_usage.to_frame().reset_index().rename(columns={'index': 'productive', 'productive': 'count'})
        productive_usage['color'] = ['#125ca4' if i==True else 'firebrick' for i in productive_usage['productive']]
        return productive_usage

    def _non_productive_distribution_plot(self):
        """
        Generate a non-productive distribution plot.
        """
        data = self._non_productive_distribution()
        chart = alt.Chart(data).mark_bar().encode(
            x=alt.X('non_productive', axis=alt.Axis(labelAngle=0), title='None productive calls'),
            y=alt.Y('count', axis=alt.Axis(format='%'), title='Percentage of reads out of non productive'),
            color=alt.Color('non_productive').scale(scheme="reds").legend(None)
        )
        return chart
    
    def _non_productive_distribution(self):
        """
        Compute non-productive distribution data.
        """
        non_productive = self.data.loc[self.data['productive'] == False].copy()
        non_productive['non_productive'] = ''
        non_productive.loc[non_productive['stop_codon'] == True, 'non_productive'] = 'stop codon present'
        non_productive.loc[non_productive['vj_in_frame'] == True, 'non_productive'] = 'VJ not in frame'
        non_productive.loc[(non_productive['vj_in_frame'] == True) & (non_productive['stop_codon'] == True), 'non_productive'] = 'stop codon present & VJ not in frame'
        counts = non_productive['non_productive'].value_counts()
        total = counts.sum()
        non_productive_usage = counts / total
        non_productive_usage = non_productive_usage.to_frame().reset_index().rename(columns={'index': 'non_productive', 'non_productive': 'count'})
        return non_productive_usage
    
    def _productive_distribution_plot_group(self):
        """
        Generate productive distribution plot group.
        """
        if self.non_productive:
            gr_plot = dp.Group(
                dp.Plot(self._productive_distribution_plot()),
                dp.Plot(self._non_productive_distribution_plot()),
                columns=2
            )
        else:
            gr_plot = dp.Plot(self._productive_distribution_plot())
        return gr_plot
    
    def _cdr3_length_distribution_plot(self):
        """
        Generate CDR3 length distribution plot.
        """
        data = self._cdr3_length_distribution()
        chart = alt.Chart(data).mark_bar().encode(x='CDR3 length', y='count')
        return chart
    
    def _cdr3_length_distribution(self):
        """
        Compute CDR3 length distribution data.
        """
        cdr3_length = self.data['junction_sequence_end'] - self.data['junction_sequence_start'] - 6
        cdr3_length = cdr3_length.value_counts()
        cdr3_length = cdr3_length.to_frame().reset_index().rename(columns={'index': 'CDR3 length', 0: 'count'})
        return cdr3_length
    
    def _v_start_distribution_plot(self):
        """
        Generate V sequence start distribution plot.
        """
        data = self._v_start_distribution()
        chart = alt.Chart(data).mark_line(interpolate="step-after").encode(x='V sequence start', y='ecdf')
        return chart
    
    def _v_start_distribution(self):
        """
        Compute V sequence start distribution data.
        """
        v_start = self.data['v_sequence_start']
        v_start = v_start.value_counts(normalize=True).sort_index().cumsum()
        v_start = v_start.to_frame().reset_index().rename(columns={'index': 'V sequence start', 'v_sequence_start': 'ecdf'})
        return v_start
    
    def _general_statistics_group(self):
        """
        Generate general statistics plot group.
        """
        if self.corrupted:
            gr_plot = dp.Group(
                dp.Plot(self._cdr3_length_distribution_plot()),
                dp.Plot(self._v_start_distribution_plot()),
                columns=2
            )
        else:
            gr_plot = dp.Plot(self._cdr3_length_distribution_plot())
        return gr_plot
    
    def _save_report(self):
        """
        Save the report using Datapane.
        """
        dp.save_report(self.layout, self.path, open=self.open)
