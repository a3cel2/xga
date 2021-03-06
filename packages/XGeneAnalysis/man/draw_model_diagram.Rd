% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nn_functions.R
\name{draw_model_diagram}
\alias{draw_model_diagram}
\title{Draw a diagram for neural network weights}
\usage{
draw_model_diagram(named_model_weights, vacuolar_transporters = c("YCF1",
  "NFT1", "YBT1"), membrane_transporters = c("PDR5", "YOR1", "AUS1", "SNQ2"),
  transporter_colors = list(PDR5 = rgb(41, 145, 206, maxColorValue = 255),
  SNQ2 = rgb(180, 231, 172, maxColorValue = 255), YOR1 = rgb(234, 151, 73,
  maxColorValue = 255), YBT1 = rgb(255, 255, 190, maxColorValue = 255), YCF1 =
  rgb(167, 0, 6, maxColorValue = 255), YOL075C = rgb(151, 165, 198,
  maxColorValue = 255), BPT1 = rgb(193, 130, 209, maxColorValue = 255), PDR12 =
  rgb(151, 165, 198, maxColorValue = 255), PDR18 = rgb(151, 165, 198,
  maxColorValue = 255), PDR15 = rgb(151, 165, 198, maxColorValue = 255), PDR11 =
  rgb(151,      165, 198, maxColorValue = 255), AUS1 = rgb(151, 165, 198,
  maxColorValue = 255), PDR10 = rgb(151, 165, 198, maxColorValue = 255), ADP1 =
  rgb(151, 165, 198, maxColorValue = 255), VMR1 = rgb(193, 130, 209,
  maxColorValue = 255), NFT1 = rgb(193, 130, 209, maxColorValue = 255)),
  drug_visual_list = list(benomyl = c(rgb(209, 66, 40, maxColorValue = 255),
  4), beauvericin = c(rgb(246, 140, 62, maxColorValue = 255), 4), cycloheximide
  = c(rgb(254, 235, 158, maxColorValue = 255), 4), methotrexate = c(rgb(127,
  116, 175, maxColorValue = 255), 3), fluconazole = c(rgb(7, 55, 86,
  maxColorValue = 255), 5), ketoconazole = c(rgb(5, 110, 174, maxColorValue =
  255), 5), itraconazole = c(rgb(126, 186, 216, maxColorValue = 255), 5),
  miconazole = c(rgb(209, 210, 230, maxColorValue = 255), 5), mitoxantrone =
  c(rgb(27,      121, 61, maxColorValue = 255), 3), camptothecin = c(rgb(166,
  213, 157, maxColorValue = 255), 3), tamoxifen = c(rgb(232, 212, 231,
  maxColorValue = 255), 3), bisantrene = c(rgb(222, 119, 172, maxColorValue =
  255), 3), cisplatin = c(rgb(145, 53, 146, maxColorValue = 255), 3), colchicine
  = c(rgb(166, 85, 40, maxColorValue = 255), 4), imatinib = c(rgb(48, 21, 78,
  maxColorValue = 255), 3), valinomycin = c(rgb(102, 8, 33, maxColorValue =
  255), 4)), big_transporters = c("PDR5", "YOR1", "SNQ2", "YBT1", "YCF1"),
  angle_positions_membrane = c(90, 170, 210, 370) * (1/180) * pi,
  angle_positions_vacuole = c(-20, 90, 200) * (1/180) * pi,
  transporter_radius = 0.15, transporter_radius_small = 0.07,
  text_size_constant = 15, max_compound_dist_from_arrow = 0.12,
  efflux_arrow_length_scale_factor = 0.03,
  efflux_arrow_width_scale_factor = 0.025,
  efflux_line_width_scale_factor = 0.33, maximum_compound_size = 0.045,
  compound_size_constant = 0.02, influence_line_width_scale_factor = 20,
  influence_line_width_max_size = Inf,
  influence_arrowhead_scale_factor = 0.05, influence_line_curvature = -0.12,
  compound_outline_scale_factor = 30, min_compound_outline_width = 0.5,
  legend_position = c(-1.25, 1.3), legend_size_left = 0.47,
  legend_size_right = 0.47, legend_text_staggering = 0.1,
  i_weights_in_legend = c(-0.04, -0.08, -0.16, -0.32, -0.64),
  e_weights_in_legend = c(0.5, 1, 2, 4, 8, 16),
  legend_influence_arrow_length = 0.3, legend_column_spacing = 0.55)
}
\arguments{
\item{named_model_weights}{a named list representing weights of the neural
network; output of name_model_weights}

\item{vacuolar_transporters}{a vector of transporters to be drawn in the vacuole}

\item{membrane_transporters}{a vector of transporters to be drawn in the membrane}

\item{transporter_colors}{a list with transporters as names and colors as values.
defines the drawn colors}

\item{drug_visual_list}{a list with drugs as names and a vector of c(color, number of sides) as values.
defines the drawn colors and shapes of the 'drugs' to be effluxed out}

\item{big_transporters}{a vector of transporters to be drawn at the 'large' size -
defaults to the five frequently-associated transporters}

\item{angle_positions_membrane}{a vector of radians, one for each transporter in membrane_transporters.
defines the angles to plot transporters on the outer membrane}

\item{angle_positions_vacuole}{a vector of radians, one for each transporter in vacuolar_transporters.
defines the angles to plot transporters on the vacuole}

\item{transporter_radius}{transporter radius (arbitrary units,roughly fraction of plot width) for
transporters in big_transporters}

\item{transporter_radius_small}{transporter radius (arbitrary units,roughly fraction of plot width) for
transporters not in big_transporters}

\item{text_size_constant}{arbitrary constant to specify text size}

\item{max_compound_dist_from_arrow}{arbitrary constant defining minimum distance from compound
drawing to the transporter drawing}

\item{efflux_arrow_length_scale_factor}{arbitrary constant defining the length of arrows representing
efflux weights (legend updated automatically)}

\item{efflux_arrow_width_scale_factor}{arbitrary constant defining the width of arrows representing
efflux weights (legend updated automatically)}

\item{efflux_line_width_scale_factor}{arbitrary constant defining the width of lines representing
efflux weights (legend updated automatically)}

\item{maximum_compound_size}{arbitrary constant defining the largest to draw compounds given
their efflux weights (legend update automatically)}

\item{compound_size_constant}{arbitrary constant defining the size of compounds given their
efflux weights (legend updated automatically)}

\item{influence_line_width_scale_factor}{arbitrary constant defining the width of
influence weight lines given their values}

\item{influence_line_width_max_size}{arbitrary constant defining the maximum width of
influence weight lines given their values}

\item{influence_arrowhead_scale_factor}{arbitrary constant defining the size of
influence weight arrows given their values}

\item{influence_line_curvature}{arbitrary constant defining the curvature
of influence weight lines}

\item{compound_outline_scale_factor}{arbitrary constant defining the outline line width of drawn
compounds, given their efflux weights (legend updated automatically)}

\item{min_compound_outline_width}{arbitrary constant defining the minimum outline line width of drawn
compounds, given their efflux weights (legend updated automatically)}

\item{legend_position}{x and y co-ordinates in the legend (plot extends from -1 to 1 in both x and y)}

\item{legend_size_left}{arbitrary constant defining the size of the left part of the legend ('I weights')}

\item{legend_size_right}{arbitrary constant defining the size of the right part of the legend ('E weights')}

\item{legend_text_staggering}{arbitrary constant specifying the gap to leave between the legend text and
the lines in the drawing}

\item{i_weights_in_legend}{a vector of I weights to display in the legend}

\item{e_weights_in_legend}{a vector of E weights to display in the legend}

\item{legend_influence_arrow_length}{arbitrary constant specifying length of influence arrows, given their weight}

\item{legend_column_spacing}{arbitrary constant specifying space between the left and right side of the legend}
}
\value{
NULL; just draws a plot
}
\description{
Draw a diagram for neural network weights
}
