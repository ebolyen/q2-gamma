const path = require('path');
const HtmlWebpackPlugin = require('html-webpack-plugin');

module.exports = {
  entry: './app/index.js',
  output: {
    filename: 'bundle.js',
    path: path.resolve(__dirname, '../dist')
  },
  mode: 'development',
  module: {
    rules: [
      {
        test: /\.rs$/,
        use: [{
          loader: 'wasm-loader'
        }, {
          loader: 'rust-native-wasm-loader',
          options: {
            release: true
          }
        }]
      }
    ]
  },
  plugins: [
    new HtmlWebpackPlugin()
  ]
}

