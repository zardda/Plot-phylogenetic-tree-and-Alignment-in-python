from pymsaviz import MsaViz
import matplotlib.pyplot as plt
from pymsaviz.config import COLOR_SCHEMES, AxesType
from matplotlib.gridspec import GridSpec
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from Bio import Phylo
from Bio.AlignIO import MultipleSeqAlignment as MSA
from Bio.SeqRecord import SeqRecord

class my_ax(MsaViz):

    def __init__(self,msa,ax,tree,ymin,ymax,color_scheme="Clustal",show_label=False):
        super().__init__(msa,color_scheme="Clustal",show_label=False)
        self.ax = ax
        self.tree = tree
        self.ymin = ymin
        self.ymax = ymax

    def get_ax(self) -> Axes:
        """Plot figure
        Parameters
        ----------
        dpi : int, optional
            Figure DPI
        Returns
        -------
        fig : Figure
            Figure
        """
        # Setup plot figure configs
        ax_type2y_size = {
            AxesType.MSA: self.msa_count * self._y_unit_size,
            AxesType.SPACE: self._y_unit_size * 1.5,
            AxesType.CONSENSUS: self._y_unit_size * self._consensus_size,
            AxesType.WRAP_SPACE: self._y_unit_size * self._wrap_space_size,
        }

        plot_ax_types = []
        for wrap_idx in range(self.wrap_num + 1):
            plot_ax_types.append(AxesType.MSA)
            if self._show_consensus:
                plot_ax_types.append(AxesType.SPACE)
                plot_ax_types.append(AxesType.CONSENSUS)
            if wrap_idx != self.wrap_num:
                plot_ax_types.append(AxesType.WRAP_SPACE)

        y_size_list = [ax_type2y_size[t] for t in plot_ax_types]
        # figsize = (self._wrap_length * self._x_unit_size, sum(y_size_list))
        # fig: Figure = plt.figure(figsize=figsize, tight_layout=True)
        # gs = GridSpec(nrows=len(plot_ax_types), ncols=1, height_ratios=y_size_list)
        # gs.update(left=0, right=1, bottom=0, top=1, hspace=0, wspace=0)

        # Plot figure
        wrap_cnt = 0
        for idx, plot_ax_type in enumerate(plot_ax_types):
            # ax: Axes = fig.add_subplot(gs[idx])
            # print(gs[idx])
            # if not isinstance(ax, Axes):
            #     raise TypeError("Error: Not matplotlib Axes class instance.")

            start = self._start + self._wrap_length * wrap_cnt
            end = self._start + self._wrap_length * (wrap_cnt + 1)
            end = self._end if end > self._end else end

            if plot_ax_type == AxesType.MSA:
                self.my_plot_msa(start, end)
            # elif plot_ax_type == AxesType.CONSENSUS:
            #     self._plot_consensus(ax, start, end)
            # elif plot_ax_type == AxesType.SPACE:
            #     ax.axis("on")
            # elif plot_ax_type == AxesType.WRAP_SPACE:
            #     ax.axis("on")
            #     wrap_cnt += 1
            # else:
            #     raise NotImplementedError(f"{plot_ax_type=} is invalid.")

        return self.ax
    
    def my_plot_msa(
        self, start, end
    ) -> None:
        """Plot MSA

        Parameters
        ----------
        ax : Axes
            Matplotlib axes to be plotted
        start : int | None, optional
            Start position. If None, `0` is set.
        end : int | None, optional
            End position. If None, `alignment_length` is set.
        """
        # Set xlim, ylim
        ax = self.ax
        start = 0 if start is None else start
        end = self.alignment_length if end is None else end
        ax.set_xlim(start, start + self._wrap_length)
        
        # reset the ylim along the tree
        ax.set_ylim(self.ymin, self.ymax)

        # Set spines & tick params (Only show bottom ticklables)
        for pos in ("left", "right", "top", "bottom"):
            ax.spines[pos].set_visible(False)
        ax.tick_params(left=False, labelleft=False)

        # Plot alignment position every 10 chars on xticks
        ticks_interval = self._ticks_interval

        # Edit the code
        # if ticks_interval is None:
            # ax.tick_params(bottom=False, labelbottom=False)
        # else:
        #     tick_ranges = range(start + 1, end + 1)
        #     xticklabels = list(filter(lambda n: n % ticks_interval == 0, tick_ranges))
        #     xticks = [n - 0.5 for n in xticklabels]
        #     ax.set_xticks(xticks, xticklabels, size=8)  # type: ignore
        ax.tick_params(bottom=False, labelbottom=False)

        plot_patches = []
        for cnt in range(self.msa_count):
            # print (cnt,y_i)
            msa_seq = self.get_sorted_seq_list()[cnt]
            y_center = self.msa_count - cnt
            y_lower = y_center - 0.5
            # Plot label text
            if self._show_label:
                if self._label_type == "id":
                    label = self.get_sorted_seq_id()[cnt]
                # elif self._label_type == "description":
                #     label = self.desc_list[cnt]
                else:
                    err_msg = f"{self._label_type=} is invalid (`id`|`description`)"
                    raise ValueError(err_msg)
                ax.text(start - 1, y_center, label, ha="right", va="center", size=10)
            # Plot count text
            if self._show_count:
                scale = end - self._start - msa_seq[self._start : end].count("-")
                ax.text(end + 1, y_center, scale, ha="left", va="center", size=10)
            for x_left in range(start, end):
                # Add colored rectangle patch
                seq_char = msa_seq[x_left]
                rect_prop = dict(
                    xy=(x_left, y_lower), width=1, height=1, color="none", lw=0
                )
                highlight_positions = self._highlight_positions
                if highlight_positions is None or x_left in highlight_positions:
                    color = self.color_scheme.get(seq_char, "#FFFFFF")
                    if self._color_scheme_name == "Identity":
                        color = self._get_identity_color(seq_char, x_left)
                    if self._custom_color_func is not None:
                        custom_color = self._custom_color_func(
                            cnt, x_left, seq_char, self.msa
                        )
                        color = color if custom_color is None else custom_color
                    rect_prop.update(**dict(color=color, lw=0, fill=True))
                if self._show_grid:
                    rect_prop.update(**dict(ec=self._grid_color, lw=0.5))
                plot_patches.append(Rectangle(**rect_prop))

                # Plot seq char text
                x_center = x_left + 0.5
                if self._show_seq_char:
                    ax.text(
                        x_center, y_center, seq_char, ha="center", va="center", size=10
                    )
                # Plot marker
                if cnt == 0 and x_left in self._pos2marker_kws:
                    marker_kws = self._pos2marker_kws[x_left]
                    ax.plot(x_center, y_center + 1, **marker_kws)
                # Plot text annotation
                if cnt == 0 and x_left in self._pos2text_kws:
                    text_kws = self._pos2text_kws[x_left]
                    ax.text(**text_kws)

        # Plot colored rectangle patch collection (Use collection for speedup)
        collection = PatchCollection(plot_patches, match_original=True, clip_on=False)
        ax.add_collection(collection)
        self.ax = ax
    
    def order_msa_id(self):
        uid2id = {}
        for idx, rec in enumerate(self._msa):
            uid = str(rec.id)
            uid2id[uid] = rec.id
            rec.id = uid
        uid2seq = {rec.id: rec.seq for rec in self._msa}
        tree = Phylo.read(self.tree,"newick")
        tree.as_phyloxml()
        tree.rooted = True
        sorted_msa = MSA([])
        for leaf in tree.get_terminals():
            uid = str(leaf.name)
            id, seq = uid2id[uid], uid2seq[uid]
            sorted_msa.append(SeqRecord(seq, id=id))
        return sorted_msa
    
    def get_sorted_seq_list(self):
        return [str(rec.seq) for rec in self.order_msa_id()]

    def get_sorted_seq_id(self):
        return [str(rec.id) for rec in self.order_msa_id()]
