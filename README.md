    gcloud builds submit --tag gcr.io/dash-plotly-testing/dash-plotly-testing_2 --project=dash-plotly-testing
    gcloud run deploy --image gcr.io/dash-plotly-testing/dash-plotly-testing_2 --platform managed --project=dash-plotly-testing --allow-unauthenticated
