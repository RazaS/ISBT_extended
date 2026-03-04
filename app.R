library(shiny)
library(DT)
library(htmltools)

`%||%` <- function(x, y) {
  if (is.null(x)) {
    return(y)
  }
  x
}

read_table <- function(path) {
  df <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = character(0),
    comment.char = ""
  )

  df[] <- lapply(df, function(col) {
    if (is.character(col)) {
      col <- gsub("\r\n", "\n", col, fixed = TRUE)
      col <- gsub("\r", "\n", col, fixed = TRUE)
      col <- gsub("_x000D_", "", col, fixed = TRUE)
    }
    col
  })

  df
}

pick_cols <- function(df, cols) {
  keep <- cols[cols %in% names(df)]
  if (length(keep) == 0) {
    return(df)
  }
  df[, keep, drop = FALSE]
}

apply_global_filter <- function(df, query) {
  q <- trimws(query %||% "")
  if (!nzchar(q) || nrow(df) == 0) {
    return(df)
  }

  q_low <- tolower(q)
  keep <- apply(df, 1, function(row_vals) {
    row_text <- tolower(paste(row_vals, collapse = " || "))
    grepl(q_low, row_text, fixed = TRUE)
  })

  df[keep, , drop = FALSE]
}

to_multiline_html <- function(x) {
  escaped <- htmlEscape(as.character(x %||% ""))
  gsub("\n", "<br/>", escaped, fixed = TRUE)
}

truncate_with_tooltip <- function(x, max_chars = 10L) {
  vapply(x, function(val) {
    if (is.null(val) || is.na(val)) {
      return("")
    }
    txt <- as.character(val)
    if (nchar(txt, type = "chars", allowNA = FALSE, keepNA = FALSE) <= max_chars) {
      return(htmlEscape(txt))
    }
    safe <- htmlEscape(txt)
    paste0("<span title=\"", safe, "\">...</span>")
  }, character(1))
}

with_row_checkboxes <- function(df) {
  n <- nrow(df)
  check_col <- if (n > 0) {
    sprintf("<input type='checkbox' class='row-check' value='%d'/>", seq_len(n))
  } else {
    character(0)
  }
  data.frame(` ` = check_col, .row_index = seq_len(n), df, check.names = FALSE, stringsAsFactors = FALSE)
}

checkbox_callback <- function(input_id, extra_js = "") {
  JS(sprintf(
    "var tbl = table;
     function sendCheckedRows() {
       var vals = [];
       tbl.$('input.row-check:checked').each(function() {
         var v = parseInt(this.value, 10);
         if (!isNaN(v)) vals.push(v);
       });
       Shiny.setInputValue('%s', vals, {priority: 'event'});
     }
     tbl.on('change', 'input.row-check', function(e) {
       e.stopPropagation();
       sendCheckedRows();
     });
     %s
     sendCheckedRows();",
    input_id,
    extra_js
  ))
}

gene_table <- read_table("Gene_table.csv")
allele_table <- read_table("Allele_table.csv")
variant_table <- read_table("Variant_table.csv")
exon_table <- read_table("Exon_table.csv")
bridge_table <- read_table("Bridge_table.csv")

allele_table$Allele_id <- as.character(allele_table$Allele_id)
variant_table$Variant_id <- as.character(variant_table$Variant_id)
bridge_table$Allele_id <- as.character(bridge_table$Allele_id)
bridge_table$Variant_id <- as.character(bridge_table$Variant_id)

bridge_by_allele <- split(bridge_table$Variant_id, bridge_table$Allele_id)
bridge_by_variant <- split(bridge_table$Allele_id, bridge_table$Variant_id)

allele_table$Variants <- vapply(allele_table$Allele_id, function(aid) {
  linked <- bridge_by_allele[[aid]]
  linked <- unique(linked[!is.na(linked) & nzchar(linked)])
  if (length(linked) == 0) {
    ""
  } else {
    paste(linked, collapse = "\n")
  }
}, character(1))

all_genes <- sort(unique(c(
  as.character(gene_table$Gene),
  as.character(allele_table$Gene),
  as.character(variant_table$Gene),
  as.character(exon_table$Gene)
)))
all_genes <- all_genes[nzchar(all_genes)]

allele_view_cols <- names(allele_table)
variant_view_cols <- names(variant_table)
exon_view_cols <- names(exon_table)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      html, body {
        height: 100%;
      }
      .container-fluid {
        padding-top: 12px;
        padding-bottom: 12px;
      }
      .app-root {
        display: flex;
        gap: 12px;
        min-height: calc(100vh - 40px);
      }
      .panel-a {
        flex: 0 0 80%;
        display: flex;
        flex-direction: column;
        gap: 12px;
        min-width: 0;
      }
      .panel-b {
        flex: 0 0 20%;
        display: flex;
        flex-direction: column;
        gap: 12px;
        min-width: 0;
      }
      .section-box {
        border: 1px solid #c8c8c8;
        border-radius: 8px;
        background: #ffffff;
        padding: 10px;
        min-width: 0;
      }
      .a1-box {
        min-height: 78px;
        display: flex;
        align-items: center;
      }
      .a2-box {
        flex: 1 1 auto;
        display: flex;
        flex-direction: column;
        min-height: 0;
        min-width: 0;
        overflow: hidden;
      }
      .mode-row {
        display: flex;
        gap: 8px;
      }
      .toggle-btn {
        width: 100%;
        border: 1px solid #b8b8b8;
        background: #f7f7f7;
      }
      .active-btn {
        background: #dbeaf8 !important;
        border-color: #6a9ecf !important;
        font-weight: 600;
      }
      .a2-allele-layout {
        display: flex;
        flex: 1 1 auto;
        max-height: 100%;
        gap: 12px;
        min-height: 0;
        min-width: 0;
        overflow: hidden;
      }
      .a2-allele-a {
        flex: 0 0 24%;
        display: grid;
        grid-template-rows: auto minmax(0, 1fr);
        max-height: 60vh;
        min-height: 0;
        min-width: 0;
        overflow: hidden;
        align-content: stretch;
      }
      .a2-allele-b {
        flex: 1 1 76%;
        display: flex;
        flex-direction: column;
        min-height: 0;
        min-width: 0;
      }
      .gene-search-wrap {
        min-height: 0;
      }
      .gene-search-wrap .form-group {
        margin-bottom: 8px;
      }
      .gene-list-wrap {
        display: flex;
        flex-direction: column;
        flex: 1 1 52vh;
        height: 52vh;
        max-height: 52vh;
        min-height: 0;
        box-sizing: border-box;
        overflow-y: auto;
        overflow-x: hidden;
        border: 1px solid #e2e2e2;
        border-radius: 6px;
        padding: 6px;
        background: #fbfbfb;
      }
      .gene-list-wrap .dataTables_wrapper,
      .gene-list-wrap .dataTables_scroll,
      .gene-list-wrap .dataTables_scrollBody {
        height: 100% !important;
        max-height: 52vh !important;
      }
      .gene-list-wrap .dataTables_scrollHead {
        display: none;
      }
      .gene-list-wrap table.dataTable tbody tr {
        cursor: pointer;
      }
      .row-check-cell {
        width: 34px;
        text-align: center;
      }
      input.row-check {
        cursor: pointer;
      }
      .detail-controls {
        flex: 0 0 auto;
      }
      .detail-controls-grid {
        display: grid;
        grid-template-columns: minmax(220px, 1fr) 160px 160px;
        gap: 8px;
        align-items: end;
      }
      .detail-controls .form-group {
        margin-bottom: 8px;
      }
      .detail-table-wrap {
        flex: 1 1 auto;
        min-height: 0;
        min-width: 0;
        overflow: hidden;
      }
      .dt-shell {
        width: 100%;
        max-width: 100%;
        min-width: 0;
        height: 100%;
        overflow-x: auto;
        overflow-y: hidden;
      }
      .dt-shell .dataTables_wrapper {
        width: 100%;
        max-width: 100%;
      }
      .dt-shell .dataTables_scroll {
        overflow-x: auto !important;
      }
      .dt-shell table.dataTable {
        max-width: none !important;
      }
      .b1-box {
        flex: 0 0 38%;
        overflow-y: auto;
      }
      .b2-box {
        flex: 1 1 62%;
        overflow-y: auto;
      }
      .last-feedback {
        margin-top: 10px;
        border-top: 1px solid #e6e6e6;
        padding-top: 8px;
        font-size: 0.92em;
      }
      .download-bar {
        margin-top: 10px;
      }
      .download-note {
        margin-top: 6px;
        font-size: 0.92em;
        color: #555;
      }
      .download-actions {
        display: flex;
        gap: 8px;
        align-items: center;
      }
      @media (max-width: 1100px) {
        .app-root {
          flex-direction: column;
          min-height: auto;
        }
        .panel-a, .panel-b {
          flex: 1 1 auto;
        }
        .a2-allele-layout {
          flex-direction: column;
        }
        .a2-allele-a, .a2-allele-b {
          flex: 1 1 auto;
        }
        .detail-controls-grid {
          grid-template-columns: 1fr;
        }
      }
    ")),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('setActiveButtons', function(msg) {
        if (!msg || !msg.groupClass) return;
        $('.' + msg.groupClass).removeClass('active-btn');
        if (msg.activeId) {
          $('#' + msg.activeId).addClass('active-btn');
        }
      });
      Shiny.addCustomMessageHandler('clearAllRowChecks', function(msg) {
        $('table.dataTable input.row-check:checked').each(function() {
          $(this).prop('checked', false).trigger('change');
        });
      });
      Shiny.addCustomMessageHandler('clearRowChecksInOutput', function(msg) {
        if (!msg || !msg.outputId) return;
        var root = $('#' + msg.outputId);
        if (!root.length) return;
        root.find('table.dataTable input.row-check:checked').each(function() {
          $(this).prop('checked', false).trigger('change');
        });
      });
      $(document).on('shiny:connected', function() {
        $('#feedback_state').prop('readonly', true);
      });
    "))
  ),
    div(
      class = "app-root",
    div(
      class = "panel-a",
      div(
        class = "section-box a1-box",
        div(
          class = "mode-row",
          actionButton("mode_gene", "Gene Table", class = "toggle-btn mode-btn active-btn"),
          actionButton("mode_allele", "Allele Table", class = "toggle-btn mode-btn"),
          actionButton("mode_exon", "Exon Table", class = "toggle-btn mode-btn")
        )
      ),
      uiOutput("panel_a2")
    ),
    div(
      class = "panel-b",
      div(
        class = "section-box download-bar",
        div(
          class = "download-actions",
          downloadButton("download_checked_csv", "Download CSV"),
          actionButton("clear_checked_rows", "Clear All", class = "btn-default")
        ),
        div(class = "download-note", textOutput("download_context"))
      ),
      div(
        class = "section-box b1-box",
        h4("Latest updates"),
        uiOutput("latest_updates")
      ),
      div(
        class = "section-box b2-box",
        h4("Feedback form"),
        textInput("feedback_state", "State", value = ""),
        textInput("feedback_email", "Email", value = ""),
        textAreaInput("feedback_comment", "Comment", value = "", rows = 5),
        actionButton("feedback_submit", "Submit", class = "btn-primary"),
        uiOutput("feedback_status")
      )
    )
  )
)

server <- function(input, output, session) {
  state <- reactiveValues(
    A2_Mode = "GT",
    Selected_Gene = if (length(all_genes) > 0) all_genes[[1]] else "",
    Gene_Search_Query = "",
    Selected_Granularity = "Allele",
    Table_Search_Query = ""
  )

  rv <- reactiveValues(
    feedback = list(),
    modal_variant_df = data.frame(),
    modal_allele_df = data.frame(),
    active_table = "gene_table",
    checked_gene = integer(0),
    checked_exon = integer(0),
    checked_detail = integer(0),
    checked_modal_variant = integer(0),
    checked_modal_allele = integer(0)
  )

  observe({
    session$sendCustomMessage(
      "setActiveButtons",
      list(
        groupClass = "mode-btn",
        activeId = if (identical(state$A2_Mode, "GT")) {
          "mode_gene"
        } else if (identical(state$A2_Mode, "AT")) {
          "mode_allele"
        } else {
          "mode_exon"
        }
      )
    )
    session$sendCustomMessage(
      "setActiveButtons",
      list(
        groupClass = "gran-btn",
        activeId = if (identical(state$Selected_Granularity, "Allele")) "gran_allele" else "gran_variant"
      )
    )
  })

  observe({
    if (!nzchar(state$Selected_Gene) || !(state$Selected_Gene %in% all_genes)) {
      state$Selected_Gene <- if (length(all_genes) > 0) all_genes[[1]] else ""
    }
  })

  observeEvent(input$mode_gene, {
    state$A2_Mode <- "GT"
    rv$active_table <- "gene_table"
  }, ignoreInit = TRUE)

  observeEvent(input$mode_allele, {
    state$A2_Mode <- "AT"
    rv$active_table <- "detail_table"
  }, ignoreInit = TRUE)

  observeEvent(input$mode_exon, {
    state$A2_Mode <- "ET"
    rv$active_table <- "exon_table_main"
  }, ignoreInit = TRUE)

  observeEvent(input$gran_allele, {
    state$Selected_Granularity <- "Allele"
  }, ignoreInit = TRUE)

  observeEvent(input$gran_variant, {
    state$Selected_Granularity <- "Variant"
  }, ignoreInit = TRUE)

  observeEvent(input$gene_search, {
    state$Gene_Search_Query <- input$gene_search %||% ""
  }, ignoreInit = FALSE)

  observeEvent(input$table_search, {
    state$Table_Search_Query <- input$table_search %||% ""
  }, ignoreInit = FALSE)

  observeEvent(list(state$Selected_Gene, state$Selected_Granularity, state$Table_Search_Query, state$A2_Mode), {
    rv$checked_detail <- integer(0)
  }, ignoreInit = TRUE)

  visible_state <- reactive({
    mode_label <- if (identical(state$A2_Mode, "GT")) {
      "Gene"
    } else if (identical(state$A2_Mode, "AT")) {
      "Allele"
    } else {
      "Exon"
    }
    gene_label <- if (identical(state$A2_Mode, "AT")) state$Selected_Gene else "N/A"
    view_label <- if (identical(state$A2_Mode, "AT")) {
      paste0(state$Selected_Granularity, "_Table")
    } else if (identical(state$A2_Mode, "ET")) {
      "Exon_Table"
    } else {
      "N/A"
    }

    paste0(
      "Mode: ", mode_label,
      " | Gene: ", gene_label,
      " | B2 View: ", view_label,
      " | GeneSearch: ", state$Gene_Search_Query,
      " | TableSearch: ", state$Table_Search_Query
    )
  })

  observe({
    updateTextInput(session, "feedback_state", value = visible_state())
  })

  filtered_genes <- reactive({
    q <- trimws(state$Gene_Search_Query)
    if (!nzchar(q)) {
      return(all_genes)
    }
    all_genes[grepl(tolower(q), tolower(all_genes), fixed = TRUE)]
  })

  observe({
    genes <- filtered_genes()
    if (length(genes) > 0 && !(state$Selected_Gene %in% genes)) {
      state$Selected_Gene <- genes[[1]]
    }
  })

  allele_for_selected_gene <- reactive({
    req(nzchar(state$Selected_Gene))
    allele_table[allele_table$Gene == state$Selected_Gene, , drop = FALSE]
  })

  variant_for_selected_gene <- reactive({
    req(nzchar(state$Selected_Gene))
    variant_table[variant_table$Gene == state$Selected_Gene, , drop = FALSE]
  })

  detail_view_raw <- reactive({
    req(identical(state$A2_Mode, "AT"))
    req(nzchar(state$Selected_Gene))

    if (identical(state$Selected_Granularity, "Allele")) {
      base_df <- pick_cols(allele_for_selected_gene(), allele_view_cols)
    } else {
      base_df <- pick_cols(variant_for_selected_gene(), variant_view_cols)
    }

    apply_global_filter(base_df, state$Table_Search_Query)
  })

  output$panel_a2 <- renderUI({
    if (identical(state$A2_Mode, "GT")) {
      div(
        class = "section-box a2-box",
        div(class = "dt-shell", DTOutput("gene_table", width = "100%"))
      )
    } else if (identical(state$A2_Mode, "ET")) {
      div(
        class = "section-box a2-box",
        div(class = "dt-shell", DTOutput("exon_table_main", width = "100%"))
      )
    } else {
      div(
        class = "section-box a2-box",
        div(
          class = "a2-allele-layout",
          div(
            class = "a2-allele-a",
            div(
              class = "gene-search-wrap",
              textInput(
                "gene_search",
                "Search genes",
                value = isolate(state$Gene_Search_Query),
                placeholder = "Type to filter genes"
              )
            ),
            div(
              class = "gene-list-wrap",
              DTOutput("gene_selector_table", width = "100%")
            )
          ),
          div(
            class = "a2-allele-b",
            div(
              class = "detail-controls",
              div(
                class = "detail-controls-grid",
                textInput(
                  "table_search",
                  "Search table",
                  value = isolate(state$Table_Search_Query),
                  placeholder = "Filter rows in active table"
                ),
                actionButton("gran_allele", "Allele", class = "toggle-btn gran-btn"),
                actionButton("gran_variant", "Variant", class = "toggle-btn gran-btn")
              )
            ),
            div(
              class = "detail-table-wrap",
              div(class = "dt-shell", DTOutput("detail_table", width = "100%"))
            )
          )
        )
      )
    }
  })

  output$gene_selector_table <- renderDT({
    genes <- filtered_genes()
    if (length(genes) == 0) {
      return(datatable(
        data.frame(Gene = character(0)),
        rownames = FALSE,
        options = list(dom = "t", paging = FALSE, searching = FALSE, info = FALSE)
      ))
    }

    gene_df <- data.frame(Gene = genes, stringsAsFactors = FALSE)
    selected_idx <- match(state$Selected_Gene, genes)

    datatable(
      gene_df,
      rownames = FALSE,
      selection = list(mode = "single", selected = if (is.na(selected_idx)) 1 else selected_idx),
      options = list(
        dom = "t",
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollY = "52vh",
        scrollCollapse = TRUE
      )
    )
  })

  observeEvent(input$gene_selector_table_rows_selected, {
    idx <- as.integer(input$gene_selector_table_rows_selected %||% integer(0))
    if (length(idx) != 1 || is.na(idx)) {
      return()
    }
    genes <- filtered_genes()
    if (idx >= 1 && idx <= length(genes)) {
      state$Selected_Gene <- genes[[idx]]
    }
  }, ignoreInit = TRUE)

  output$gene_table <- renderDT({
    view_df <- with_row_checkboxes(gene_table)
    row_idx_col <- match(".row_index", names(view_df)) - 1
    callback <- checkbox_callback("gene_table_checked_rows")
    datatable(
      view_df,
      rownames = FALSE,
      selection = "none",
      callback = callback,
      escape = setdiff(names(view_df), " "),
      options = list(
        paging = FALSE,
        scrollX = TRUE,
        scrollY = "68vh",
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = 0, className = "row-check-cell", orderable = FALSE),
          list(targets = row_idx_col, visible = FALSE)
        )
      )
    )
  })

  output$exon_table_main <- renderDT({
    view_df <- with_row_checkboxes(exon_table)
    row_idx_col <- match(".row_index", names(view_df)) - 1
    callback <- checkbox_callback("exon_table_main_checked_rows")
    datatable(
      view_df,
      rownames = FALSE,
      selection = "none",
      callback = callback,
      escape = setdiff(names(view_df), " "),
      options = list(
        paging = FALSE,
        scrollX = TRUE,
        scrollY = "68vh",
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = 0, className = "row-check-cell", orderable = FALSE),
          list(targets = row_idx_col, visible = FALSE)
        )
      )
    )
  })

  output$detail_table <- renderDT({
    req(identical(state$A2_Mode, "AT"))
    req(nzchar(state$Selected_Gene))

    if (identical(state$Selected_Granularity, "Allele")) {
      raw_df <- detail_view_raw()
      display_df <- raw_df
      display_df$Variants <- to_multiline_html(display_df$Variants)
      view_df <- with_row_checkboxes(display_df)

      id_col <- match("Allele_id", names(view_df)) - 1
      row_idx_col <- match(".row_index", names(view_df)) - 1
      escape_cols <- setdiff(names(view_df), c(" ", "Variants"))
      callback <- checkbox_callback(
        "detail_table_checked_rows",
        sprintf(
          "tbl.on('dblclick', 'tbody td', function(e) {
             if (this.cellIndex === 0 || $(e.target).is('input.row-check')) return;
             var row = tbl.row(this.parentNode).data();
             if (!row) return;
             Shiny.setInputValue('allele_row_dblclick', row[%d], {priority: 'event'});
           });",
          id_col
        )
      )

      datatable(
        view_df,
        rownames = FALSE,
        selection = "none",
        escape = escape_cols,
        callback = callback,
        options = list(
          searching = FALSE,
          paging = FALSE,
          scrollX = TRUE,
          scrollY = "52vh",
          scrollCollapse = TRUE,
          autoWidth = TRUE,
          columnDefs = list(
            list(targets = 0, className = "row-check-cell", orderable = FALSE),
            list(targets = row_idx_col, visible = FALSE)
          )
        )
      )
    } else {
      raw_df <- detail_view_raw()
      display_df <- raw_df
      if ("Ref_allele_curated" %in% names(display_df)) {
        display_df$Ref_allele_curated <- truncate_with_tooltip(display_df$Ref_allele_curated, 10L)
      }
      if ("Alt_allele_curated" %in% names(display_df)) {
        display_df$Alt_allele_curated <- truncate_with_tooltip(display_df$Alt_allele_curated, 10L)
      }
      view_df <- with_row_checkboxes(display_df)

      id_col <- match("Variant_id", names(view_df)) - 1
      row_idx_col <- match(".row_index", names(view_df)) - 1
      callback <- checkbox_callback(
        "detail_table_checked_rows",
        sprintf(
          "tbl.on('click', 'tbody td', function(e) {
             if (this.cellIndex === 0 || $(e.target).is('input.row-check')) return;
             var row = tbl.row(this.parentNode).data();
             if (!row) return;
             Shiny.setInputValue('variant_row_click', row[%d], {priority: 'event'});
           });",
          id_col
        )
      )

      datatable(
        view_df,
        rownames = FALSE,
        selection = "none",
        callback = callback,
        escape = setdiff(names(view_df), c(" ", "Ref_allele_curated", "Alt_allele_curated")),
        options = list(
          searching = FALSE,
          paging = FALSE,
          scrollX = TRUE,
          scrollY = "52vh",
          scrollCollapse = TRUE,
          autoWidth = TRUE,
          columnDefs = list(
            list(targets = 0, className = "row-check-cell", orderable = FALSE),
            list(targets = row_idx_col, visible = FALSE)
          )
        )
      )
    }
  })

  observeEvent(input$allele_row_dblclick, {
    allele_id <- as.character(input$allele_row_dblclick %||% "")
    req(nzchar(allele_id))

    linked_variant_ids <- unique(bridge_by_allele[[allele_id]])
    linked_variant_ids <- linked_variant_ids[!is.na(linked_variant_ids) & nzchar(linked_variant_ids)]

    modal_df <- variant_table[variant_table$Variant_id %in% linked_variant_ids, , drop = FALSE]
    if (nrow(modal_df) > 0 && length(linked_variant_ids) > 0) {
      modal_df <- modal_df[match(linked_variant_ids, modal_df$Variant_id), , drop = FALSE]
      modal_df <- modal_df[!is.na(modal_df$Variant_id), , drop = FALSE]
    }
    modal_df <- pick_cols(modal_df, variant_view_cols)
    rv$modal_variant_df <- modal_df
    rv$checked_modal_variant <- integer(0)

    allele_name <- allele_table$Allele_name[match(allele_id, allele_table$Allele_id)]
    allele_label <- if (length(allele_name) > 0 && nzchar(allele_name[[1]])) {
      allele_name[[1]]
    } else {
      allele_id
    }

    showModal(
      modalDialog(
        title = paste0("Variants associated with allele ", allele_label),
        div(class = "dt-shell", DTOutput("modal_variant_table", width = "100%")),
        easyClose = TRUE,
        size = "l",
        footer = tagList(
          downloadButton("download_modal_variant_csv", "Download CSV"),
          actionButton("clear_modal_variant_checked", "Clear All", class = "btn-default"),
          modalButton("Close")
        )
      )
    )
  }, ignoreInit = TRUE)

  observeEvent(input$variant_row_click, {
    variant_id <- as.character(input$variant_row_click %||% "")
    req(nzchar(variant_id))

    linked_allele_ids <- unique(bridge_by_variant[[variant_id]])
    linked_allele_ids <- linked_allele_ids[!is.na(linked_allele_ids) & nzchar(linked_allele_ids)]

    modal_df <- allele_table[allele_table$Allele_id %in% linked_allele_ids, , drop = FALSE]
    if (nrow(modal_df) > 0 && length(linked_allele_ids) > 0) {
      modal_df <- modal_df[match(linked_allele_ids, modal_df$Allele_id), , drop = FALSE]
      modal_df <- modal_df[!is.na(modal_df$Allele_id), , drop = FALSE]
    }
    modal_df <- pick_cols(modal_df, allele_view_cols)
    rv$modal_allele_df <- modal_df
    rv$checked_modal_allele <- integer(0)

    variant_notation <- variant_table$Nucleotide_change[match(variant_id, variant_table$Variant_id)]
    variant_label <- if (length(variant_notation) > 0 && nzchar(variant_notation[[1]])) {
      paste0(variant_id, " (", variant_notation[[1]], ")")
    } else {
      variant_id
    }

    showModal(
      modalDialog(
        title = paste0("Alleles associated with variant ", variant_label),
        div(class = "dt-shell", DTOutput("modal_allele_table", width = "100%")),
        easyClose = TRUE,
        size = "l",
        footer = tagList(
          downloadButton("download_modal_allele_csv", "Download CSV"),
          actionButton("clear_modal_allele_checked", "Clear All", class = "btn-default"),
          modalButton("Close")
        )
      )
    )
  }, ignoreInit = TRUE)

  output$modal_variant_table <- renderDT({
    modal_df <- rv$modal_variant_df
    if ("Ref_allele_curated" %in% names(modal_df)) {
      modal_df$Ref_allele_curated <- truncate_with_tooltip(modal_df$Ref_allele_curated, 10L)
    }
    if ("Alt_allele_curated" %in% names(modal_df)) {
      modal_df$Alt_allele_curated <- truncate_with_tooltip(modal_df$Alt_allele_curated, 10L)
    }
    df <- with_row_checkboxes(modal_df)
    row_idx_col <- match(".row_index", names(df)) - 1
    callback <- checkbox_callback("modal_variant_table_checked_rows")
    datatable(
      df,
      rownames = FALSE,
      selection = "none",
      callback = callback,
      escape = setdiff(names(df), c(" ", "Ref_allele_curated", "Alt_allele_curated")),
      options = list(
        paging = FALSE,
        scrollX = TRUE,
        scrollY = "45vh",
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = 0, className = "row-check-cell", orderable = FALSE),
          list(targets = row_idx_col, visible = FALSE)
        )
      )
    )
  })

  output$modal_allele_table <- renderDT({
    df <- rv$modal_allele_df
    if ("Variants" %in% names(df)) {
      df$Variants <- to_multiline_html(df$Variants)
    }
    df <- with_row_checkboxes(df)
    row_idx_col <- match(".row_index", names(df)) - 1
    callback <- checkbox_callback("modal_allele_table_checked_rows")
    datatable(
      df,
      rownames = FALSE,
      selection = "none",
      callback = callback,
      escape = setdiff(names(df), c(" ", "Variants")),
      options = list(
        paging = FALSE,
        scrollX = TRUE,
        scrollY = "45vh",
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = 0, className = "row-check-cell", orderable = FALSE),
          list(targets = row_idx_col, visible = FALSE)
        )
      )
    )
  })

  observeEvent(input$gene_table_checked_rows, {
    rv$checked_gene <- as.integer(input$gene_table_checked_rows %||% integer(0))
    rv$active_table <- "gene_table"
  }, ignoreInit = TRUE)

  observeEvent(input$exon_table_main_checked_rows, {
    rv$checked_exon <- as.integer(input$exon_table_main_checked_rows %||% integer(0))
    rv$active_table <- "exon_table_main"
  }, ignoreInit = TRUE)

  observeEvent(input$detail_table_checked_rows, {
    rv$checked_detail <- as.integer(input$detail_table_checked_rows %||% integer(0))
    rv$active_table <- "detail_table"
  }, ignoreInit = TRUE)

  observeEvent(input$modal_variant_table_checked_rows, {
    rv$checked_modal_variant <- as.integer(input$modal_variant_table_checked_rows %||% integer(0))
    rv$active_table <- "modal_variant_table"
  }, ignoreInit = TRUE)

  observeEvent(input$modal_allele_table_checked_rows, {
    rv$checked_modal_allele <- as.integer(input$modal_allele_table_checked_rows %||% integer(0))
    rv$active_table <- "modal_allele_table"
  }, ignoreInit = TRUE)

  observeEvent(input$clear_checked_rows, {
    rv$checked_gene <- integer(0)
    rv$checked_exon <- integer(0)
    rv$checked_detail <- integer(0)
    rv$checked_modal_variant <- integer(0)
    rv$checked_modal_allele <- integer(0)
    session$sendCustomMessage("clearAllRowChecks", list())
  }, ignoreInit = TRUE)

  observeEvent(input$clear_modal_variant_checked, {
    rv$checked_modal_variant <- integer(0)
    session$sendCustomMessage("clearRowChecksInOutput", list(outputId = "modal_variant_table"))
  }, ignoreInit = TRUE)

  observeEvent(input$clear_modal_allele_checked, {
    rv$checked_modal_allele <- integer(0)
    session$sendCustomMessage("clearRowChecksInOutput", list(outputId = "modal_allele_table"))
  }, ignoreInit = TRUE)

  get_active_table_label <- reactive({
    switch(
      rv$active_table,
      gene_table = "Gene Table",
      exon_table_main = "Exon Table",
      detail_table = if (identical(state$Selected_Granularity, "Allele")) "Allele Detail Table" else "Variant Detail Table",
      modal_variant_table = "Variants Popup Table",
      modal_allele_table = "Alleles Popup Table",
      "Unknown"
    )
  })

  get_active_rows <- reactive({
    table_id <- rv$active_table
    checked_rows <- switch(
      table_id,
      gene_table = rv$checked_gene,
      exon_table_main = rv$checked_exon,
      detail_table = rv$checked_detail,
      modal_variant_table = rv$checked_modal_variant,
      modal_allele_table = rv$checked_modal_allele,
      integer(0)
    )

    source_df <- switch(
      table_id,
      gene_table = gene_table,
      exon_table_main = exon_table,
      detail_table = if (identical(state$A2_Mode, "AT")) detail_view_raw() else data.frame(),
      modal_variant_table = rv$modal_variant_df,
      modal_allele_table = rv$modal_allele_df,
      data.frame()
    )

    list(
      table_id = table_id,
      label = get_active_table_label(),
      selected_rows = checked_rows,
      source_df = source_df
    )
  })

  output$download_context <- renderText({
    info <- get_active_rows()
    paste0(
      "Active table: ", info$label,
      " | Checked rows: ", length(info$selected_rows)
    )
  })

  output$download_checked_csv <- downloadHandler(
    filename = function() {
      info <- get_active_rows()
      paste0(
        gsub("[^A-Za-z0-9]+", "_", tolower(info$label)),
        "_checked_rows_",
        format(Sys.Date(), "%Y%m%d"),
        ".csv"
      )
    },
    content = function(file) {
      info <- get_active_rows()
      idx <- as.integer(info$selected_rows)
      df <- info$source_df

      if (nrow(df) == 0 || length(idx) == 0) {
        write.csv(data.frame(), file, row.names = FALSE, na = "")
        return()
      }

      idx <- unique(idx[!is.na(idx)])
      idx <- idx[idx >= 1 & idx <= nrow(df)]
      if (length(idx) == 0) {
        write.csv(data.frame(), file, row.names = FALSE, na = "")
        return()
      }

      write.csv(df[idx, , drop = FALSE], file, row.names = FALSE, na = "")
    }
  )

  output$download_modal_variant_csv <- downloadHandler(
    filename = function() {
      paste0("variants_popup_checked_rows_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      idx <- unique(as.integer(rv$checked_modal_variant))
      idx <- idx[!is.na(idx)]
      df <- rv$modal_variant_df

      if (nrow(df) == 0 || length(idx) == 0) {
        write.csv(data.frame(), file, row.names = FALSE, na = "")
        return()
      }

      idx <- idx[idx >= 1 & idx <= nrow(df)]
      if (length(idx) == 0) {
        write.csv(data.frame(), file, row.names = FALSE, na = "")
        return()
      }

      write.csv(df[idx, , drop = FALSE], file, row.names = FALSE, na = "")
    }
  )

  output$download_modal_allele_csv <- downloadHandler(
    filename = function() {
      paste0("alleles_popup_checked_rows_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      idx <- unique(as.integer(rv$checked_modal_allele))
      idx <- idx[!is.na(idx)]
      df <- rv$modal_allele_df

      if (nrow(df) == 0 || length(idx) == 0) {
        write.csv(data.frame(), file, row.names = FALSE, na = "")
        return()
      }

      idx <- idx[idx >= 1 & idx <= nrow(df)]
      if (length(idx) == 0) {
        write.csv(data.frame(), file, row.names = FALSE, na = "")
        return()
      }

      write.csv(df[idx, , drop = FALSE], file, row.names = FALSE, na = "")
    }
  )

  output$latest_updates <- renderUI({
    tags$div(
      tags$ul(
        tags$li(sprintf("Loaded Gene_table: %s rows", format(nrow(gene_table), big.mark = ","))),
        tags$li(sprintf("Loaded Allele_table: %s rows", format(nrow(allele_table), big.mark = ","))),
        tags$li(sprintf("Loaded Variant_table: %s rows", format(nrow(variant_table), big.mark = ","))),
        tags$li(sprintf("Loaded Exon_table: %s rows", format(nrow(exon_table), big.mark = ","))),
        tags$li(sprintf("Loaded Bridge_table: %s rows", format(nrow(bridge_table), big.mark = ","))),
        tags$li("Use A1 buttons to switch between Gene, Allele, and Exon exploration modes."),
        tags$li("Allele rows: double-click to inspect linked variants."),
        tags$li("Variant rows: single-click to inspect linked alleles.")
      )
    )
  })

  observeEvent(input$feedback_submit, {
    entry <- list(
      submitted_at = Sys.time(),
      state = visible_state(),
      email = input$feedback_email %||% "",
      comment = input$feedback_comment %||% ""
    )

    rv$feedback <- c(rv$feedback, list(entry))

    updateTextInput(session, "feedback_email", value = "")
    updateTextAreaInput(session, "feedback_comment", value = "")

    showNotification("Feedback submitted.", type = "message")
  }, ignoreInit = TRUE)

  output$feedback_status <- renderUI({
    if (length(rv$feedback) == 0) {
      return(tags$div(class = "last-feedback", "No submissions yet."))
    }

    last <- rv$feedback[[length(rv$feedback)]]
    tags$div(
      class = "last-feedback",
      tags$strong("Last submission"),
      tags$div(format(last$submitted_at, "%Y-%m-%d %H:%M:%S")),
      tags$div(paste0("Email: ", if (nzchar(last$email)) last$email else "(blank)")),
      tags$div(paste0("State: ", last$state)),
      tags$div(paste0("Comment length: ", nchar(last$comment)))
    )
  })
}

shinyApp(ui, server)
